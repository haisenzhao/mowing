#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"


namespace hpcg {
	
	void ToolpathGenerator::BuildeImageSpace()
	{
		//boundary_boundingbox
		double min_x = contours.bbox().xmin();
		double min_y = contours.bbox().ymin();
		double max_x = contours.bbox().xmax();
		double max_y = contours.bbox().ymax();

		//image space
		int pixel_nb = ((2 * image_space.pixel_size + max_x - min_x) / image_space.pixel_size)* ((2 * image_space.pixel_size + max_y - min_y) / image_space.pixel_size);


		int pixel_index = 0;
		for (double x = min_x - image_space.pixel_size; x < max_x + image_space.pixel_size; x += image_space.pixel_size)
		{
			std::vector<Pixel> pixels;

			for (double y = min_y - image_space.pixel_size; y < max_y + image_space.pixel_size; y += image_space.pixel_size)
			{
				pixel_index++;
				if (pixel_index % ((int)pixel_nb / 100) == 0)
				{
					std::cout << (double)pixel_index / pixel_nb << std::endl;
				}
				Pixel pixel(Vector2d(x, y));
				
				if (contours.outer_boundary().bounded_side(Point_2(x, y)) == CGAL::ON_BOUNDED_SIDE)
				{
					pixel.inside = true;

					for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
					{
						if (hole_iter->bounded_side(Point_2(x, y)) == CGAL::ON_BOUNDED_SIDE)
						{
							pixel.inside = false;
							break;
						}
					}
				}
				else
				{
					pixel.inside = false;
				}

				pixels.push_back(pixel);
			}
			image_space.pixels.push_back(pixels);
			std::vector<Pixel>().swap(pixels);
		}
	}



	void ToolpathGenerator::DrawPixel(int i, int j)
	{
		if (image_space.check(i, j))
		{
			glVertex3f(image_space.pixels[i][j].center[0] - image_space.pixel_size / 2.0, image_space.pixels[i][j].center[1] - image_space.pixel_size / 2.0, 0.0);
			glVertex3f(image_space.pixels[i][j].center[0] - image_space.pixel_size / 2.0, image_space.pixels[i][j].center[1] + image_space.pixel_size / 2.0, 0.0);
			glVertex3f(image_space.pixels[i][j].center[0] + image_space.pixel_size / 2.0, image_space.pixels[i][j].center[1] + image_space.pixel_size / 2.0, 0.0);
			glVertex3f(image_space.pixels[i][j].center[0] + image_space.pixel_size / 2.0, image_space.pixels[i][j].center[1] - image_space.pixel_size / 2.0, 0.0);
		}
	}

}