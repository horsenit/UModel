#include "UnTextureNVTT.h"

#include "Core.h"
#include "UnCore.h"
#include "UnObject.h"
#include "UnMaterial.h"
#include "UnMaterial3.h"

#include "Exporters.h"

#include "libpng/png.h"

#define TGA_SAVE_BOTTOMLEFT	1


#define TGA_ORIGIN_MASK		0x30
#define TGA_BOTLEFT			0x00
#define TGA_BOTRIGHT		0x10					// unused
#define TGA_TOPLEFT			0x20
#define TGA_TOPRIGHT		0x30					// unused

#if _MSC_VER
#pragma pack(push,1)
#endif

struct GCC_PACK tgaHdr_t
{
	byte 	id_length, colormap_type, image_type;
	uint16	colormap_index, colormap_length;
	byte	colormap_size;
	uint16	x_origin, y_origin;				// unused
	uint16	width, height;
	byte	pixel_size, attributes;
};

#if _MSC_VER
#pragma pack(pop)
#endif


bool GNoTgaCompress = false;
bool GExportDDS = false;
bool GExportPNG = false;

struct PngWriteCtx
{
	FArchive *Ar;
};

static void user_write_fn(png_structp png_ptr, png_bytep data, png_size_t length)
{
	PngWriteCtx* ctx = (PngWriteCtx*)png_get_io_ptr(png_ptr);
	ctx->Ar->Serialize(data, length);
}

static void user_flush_fn(png_structp png_ptr)
{
}

static void user_error_fn(png_structp png_ptr, png_const_charp error_msg)
{
	appError("Error writing PNG: %s", error_msg);
}

static void user_warning_fn(png_structp png_ptr, png_const_charp warning_msg)
{
	appNotify("Warning writing PNG: %s", warning_msg);
}

static void* user_malloc(png_structp /*png_ptr*/, png_size_t size)
{
	return appMalloc(size);
}

static void user_free(png_structp /*png_ptr*/, png_voidp struct_ptr)
{
	appFree(struct_ptr);
}

void WritePNG(FArchive &Ar, int width, int height, byte *pic)
{
	guard(WritePNG);

	PngWriteCtx ctx;
	ctx.Ar = &Ar;

	png_structp png_ptr = png_create_write_struct_2(PNG_LIBPNG_VER_STRING, NULL, user_error_fn, user_warning_fn, NULL, user_malloc, user_free);
	assert(png_ptr);
	png_infop info_ptr = png_create_info_struct(png_ptr);
	assert(info_ptr);
	png_set_write_fn(png_ptr, &ctx, user_write_fn, user_flush_fn);

	int size = width * height;
	int i;
	byte *src;

	// check for grayscale and alpha channel
	bool gray = true;
	bool alpha = false;
	for (i = 0, src = pic; i < size; i++, src += 4)
	{
		if (src[0] != src[1] || src[0] != src[2])
			gray = false;
		if (src[3] != 255)
			alpha = true;
		if (!gray && alpha)
			break;
	}

	png_set_compression_level(png_ptr, 4);
	png_set_IHDR(png_ptr, info_ptr, width, height, 8,
		!gray ? (alpha ? PNG_COLOR_TYPE_RGB_ALPHA : PNG_COLOR_TYPE_RGB) :
				(alpha ? PNG_COLOR_TYPE_GRAY_ALPHA : PNG_COLOR_TYPE_GRAY),
		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png_ptr, info_ptr);

	int channels = (!gray ? (alpha ? 4 : 3) : (alpha ? 2 : 1));
	int rowSize = width * channels;
	png_byte *row = (png_byte*)appMalloc(rowSize);
	for (i = 0; i < height; i++) {
		src = pic + i * (width * 4);
		byte *dst = row;
		for (int x = 0; x < width; x++) {
			if (gray) {
				*dst = *src;
				dst++;
			} else {
				*(dst+0) = *(src+0);
				*(dst+1) = *(src+1);
				*(dst+2) = *(src+2);
				dst += 3;
			}
			if (alpha) {
				*dst = *(src+3);
				dst++;
			}
			src += 4;
		}
		png_write_row(png_ptr, row);
	}
	appFree(row);

	png_write_end(png_ptr, info_ptr);
	png_destroy_write_struct(&png_ptr, NULL);
	unguard;
}

//?? place this function outside (cannot place to Core - using FArchive)

void WriteTGA(FArchive &Ar, int width, int height, byte *pic)
{
	guard(WriteTGA);

	int		i;

	byte *src;
	int size = width * height;
	// convert RGB to BGR (inplace!)
	for (i = 0, src = pic; i < size; i++, src += 4)
		Exchange(src[0], src[2]);

	// check for 24 bit image possibility
	int colorBytes = 3;
	for (i = 0, src = pic + 3; i < size; i++, src += 4)
		if (src[0] != 255)									// src initialized with offset 3
		{
			colorBytes = 4;
			break;
		}

	byte *packed = (byte*)appMalloc(width * height * colorBytes);
	byte *threshold = packed + width * height * colorBytes - 16; // threshold for "dst"

	src = pic;
	byte *dst = packed;
	int column = 0;
	byte *flag = NULL;
	bool rle = false;

	bool useCompression = true;
	if (GNoTgaCompress)
		threshold = dst;									// will break compression loop immediately

	for (i = 0; i < size; i++)
	{
		if (dst >= threshold)								// when compressed is too large, same uncompressed
		{
			useCompression = false;
			break;
		}

		byte b = *src++;
		byte g = *src++;
		byte r = *src++;
		byte a = *src++;

		if (column < width - 1 &&							// not on screen edge; NOTE: when i == size-1, col==width-1
			b == src[0] && g == src[1] && r == src[2] && a == src[3] &&	// next pixel will be the same
			!(rle && flag && *flag == 254))					// flag overflow
		{
			if (!rle || !flag)
			{
				// starting new RLE sequence
				flag = dst++;
				*flag = 128 - 1;							// will be incremented below
				*dst++ = b; *dst++ = g; *dst++ = r;			// store RGB
				if (colorBytes == 4) *dst++ = a;			// store alpha
			}
			(*flag)++;										// enqueue one more texel
			rle = true;
		}
		else
		{
			if (rle)
			{
				// previous block was RLE, and next (now: current) byte was
				// the same - enqueue it to previous block and close block
				(*flag)++;
				flag = NULL;
			}
			else
			{
				if (column == 0)							// check for screen edge
					flag = NULL;
				if (!flag)
				{
					// start new copy sequence
					flag = dst++;
					*flag = 255;
				}
				*dst++ = b; *dst++ = g; *dst++ = r;			// store RGB
				if (colorBytes == 4) *dst++ = a;			// store alpha
				(*flag)++;
				if (*flag == 127) flag = NULL;				// check for overflow
			}
			rle = false;
		}

		if (++column == width) column = 0;
	}

	// write header
	tgaHdr_t header;
	memset(&header, 0, sizeof(header));
	header.width  = width;
	header.height = height;
#if 0
	// debug: write black/white image
	header.pixel_size = 8;
	header.image_type = 3;
	fwrite(&header, 1, sizeof(header), f);
	for (i = 0; i < width * height; i++, pic += 4)
	{
		int c = (pic[0]+pic[1]+pic[2]) / 3;
		fwrite(&c, 1, 1, f);
	}
#else
	header.pixel_size = colorBytes * 8;
#if TGA_SAVE_BOTTOMLEFT
	header.attributes = TGA_BOTLEFT;
#else
	header.attributes = TGA_TOPLEFT;
#endif
	if (useCompression)
	{
		header.image_type = 10;		// RLE
		// write data
		Ar.Serialize(&header, sizeof(header));
		Ar.Serialize(packed, dst - packed);
	}
	else
	{
		header.image_type = 2;		// uncompressed
		// convert to 24 bits image, when needed
		if (colorBytes == 3)
			for (i = 0, src = dst = pic; i < size; i++)
			{
				*dst++ = *src++;
				*dst++ = *src++;
				*dst++ = *src++;
				src++;
			}
		// write data
		Ar.Serialize(&header, sizeof(header));
		Ar.Serialize(pic, size * colorBytes);
	}
#endif

	appFree(packed);

	unguard;
}

// Radiance file format

static void float2rgbe(float red, float green, float blue, byte* rgbe)
{
	float v;
	int e;
	v = red;
	if (green > v) v = green;
	if (blue > v) v = blue;
	if (v < 1e-32)
	{
		rgbe[0] = rgbe[1] = rgbe[2] = rgbe[3] = 0;
	}
	else
	{
		v = frexp(v, &e) * 256.0f / v; //?? TODO: check if frexp is slow and could be optimized
		rgbe[0] = byte(red * v);
		rgbe[1] = byte(green * v);
		rgbe[2] = byte(blue * v);
		rgbe[3] = byte(e + 128);
	}
}

static void WriteHDR(FArchive &Ar, int width, int height, byte *pic)
{
	guard(WriteHDR);

	char hdr[64];
	appSprintf(ARRAY_ARG(hdr), "#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y %d +X %d\n", height, width);

	//!! TODO: compress HDR file (seems have RLE support)
	// Convert float[w*h*4] to rgbe[w*h] inplace
	const float* floatSrc = (float*)pic;
	byte* byteDst = pic;
	for (int p = 0; p < width * height; p++, floatSrc += 4, byteDst += 4)
	{
		float2rgbe(floatSrc[0], floatSrc[1], floatSrc[2], byteDst);
	}

	Ar.Serialize(hdr, strlen(hdr));
	Ar.Serialize(pic, width * height * 4);

	unguard;
}


static void WriteDDS(const CTextureData &TexData, const char *Filename)
{
	guard(WriteDDS);

	if (!TexData.Mips.Num()) return;
	const CMipMap& Mip = TexData.Mips[0];

	unsigned fourCC = TexData.GetFourCC();

	// code from CTextureData::Decompress()
	if (!fourCC)
		appError("unknown texture format %d \n", TexData.Format);	// should not happen - IsDXT() should not pass execution here

	nv::DDSHeader header;
	header.setFourCC(fourCC & 0xFF, (fourCC >> 8) & 0xFF, (fourCC >> 16) & 0xFF, (fourCC >> 24) & 0xFF);
//	header.setPixelFormat(32, 0xFF, 0xFF << 8, 0xFF << 16, 0xFF << 24);	// bit count and per-channel masks
	//!! Note: should use setFourCC for compressed formats, and setPixelFormat for uncompressed - these functions are
	//!! incompatible. When fourcc is used, color masks are zero, and vice versa.
	header.setWidth(Mip.USize);
	header.setHeight(Mip.VSize);
//	header.setNormalFlag(TexData.Format == TPF_DXT5N || TexData.Format == TPF_3DC); -- required for decompression only
	header.setLinearSize(Mip.DataSize);

	appMakeDirectoryForFile(Filename);

	byte headerBuffer[128];							// DDS header is 128 bytes long
	memset(headerBuffer, 0, 128);
	WriteDDSHeader(headerBuffer, header);
	FArchive *Ar = new FFileWriter(Filename);
	if (!Ar->IsOpen())
	{
		appPrintf("Error creating file \"%s\" ...\n", Filename);
		delete Ar;
		return;
	}
	Ar->Serialize(headerBuffer, 128);
	Ar->Serialize(const_cast<byte*>(Mip.CompressedData), Mip.DataSize);
	delete Ar;

	unguard;
}

void ExportTexturePNGArchive(const UUnrealMaterial *Tex, FArchive &Ar, ExportPNGCommand command)
{
	guard(ExportTexturePNGArchive);

	byte *pic = NULL;
	int width, height;

	CTextureData TexData;
	if (Tex->GetTextureData(TexData))
	{
		width = TexData.Mips[0].USize;
		height = TexData.Mips[0].VSize;
		pic = TexData.Decompress();
	}

	// HDR not supported (?)
	if (!pic || PixelFormatInfo[TexData.Format].Float)
	{
		appPrintf("WARNING: texture %s has no valid mipmaps\n", Tex->Name);
		// produce 1x1-pixel tga
		// should erase file?
		width = height = 1;
		pic = new byte[4];
	}

	if (command & (ExportPNG_FlipR|ExportPNG_FlipG|ExportPNG_FlipB|ExportPNG_FlipA))
	{
		byte *end = pic + (height * width) * 4;
		for (byte *src = pic; src < end; src += 4)
		{
			if (command & ExportPNG_FlipR)
				*(src+0) = 255 - *(src+0);
			if (command & ExportPNG_FlipG)
				*(src+1) = 255 - *(src+1);
			if (command & ExportPNG_FlipB)
				*(src+2) = 255 - *(src+2);
			if (command & ExportPNG_FlipA)
				*(src+3) = 255 - *(src+3);
		}
	}

	WritePNG(Ar, width, height, pic);

	delete pic;

	Tex->ReleaseTextureData();

	unguard;
}

void ExportTexture(const UUnrealMaterial *Tex)
{
	guard(ExportTexture);

	if (GDontOverwriteFiles)
	{
		if (CheckExportFilePresence(Tex, "%s.tga", Tex->Name)) return;
		if (CheckExportFilePresence(Tex, "%s.dds", Tex->Name)) return;
		if (CheckExportFilePresence(Tex, "%s.png", Tex->Name)) return;
		if (CheckExportFilePresence(Tex, "%s.hdr", Tex->Name)) return;
	}

	byte *pic = NULL;
	int width, height;

	CTextureData TexData;
	if (Tex->GetTextureData(TexData))
	{
		if (GExportDDS && TexData.IsDXT())
		{
			WriteDDS(TexData, GetExportFileName(Tex, "%s.dds", Tex->Name));
			return;
		}

		width = TexData.Mips[0].USize;
		height = TexData.Mips[0].VSize;
		pic = TexData.Decompress();
	}

	if (!pic)
	{
		appPrintf("WARNING: texture %s has no valid mipmaps\n", Tex->Name);
		// produce 1x1-pixel tga
		// should erase file?
		width = height = 1;
		pic = new byte[4];
	}

	// For HDR textures use Radiance format
	if (PixelFormatInfo[TexData.Format].Float)
	{
		FArchive *Ar = CreateExportArchive(Tex, 0, "%s.hdr", Tex->Name);
		if (Ar)
		{
			WriteHDR(*Ar, width, height, pic);
			delete Ar;
		}

		delete pic;
		return;
	}

	if (!GExportPNG)
	{
		FArchive *Ar = CreateExportArchive(Tex, 0, "%s.tga", Tex->Name);
		if (Ar)
		{
	#if TGA_SAVE_BOTTOMLEFT
			// flip image vertically (UnrealEd for UE2 have a bug with importing TGA_TOPLEFT images,
			// it simply ignores orientation flags)
			for (int i = 0; i < height / 2; i++)
			{
				byte *p1 = pic + width * 4 * i;
				byte *p2 = pic + width * 4 * (height - i - 1);
				for (int j = 0; j < width * 4; j++)
					Exchange(p1[j], p2[j]);
			}
	#endif
			WriteTGA(*Ar, width, height, pic);
			delete Ar;
		}
	}
	else
	{
		FArchive *Ar = CreateExportArchive(Tex, 0, "%s.png", Tex->Name);
		if (Ar)
		{
			WritePNG(*Ar, width, height, pic);
			delete Ar;
		}
	}

	delete pic;

	Tex->ReleaseTextureData();

	unguard;
}
