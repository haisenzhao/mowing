#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"


namespace hpcg {
	
	bool goon = false;
	bool over = false;
	int iii_debug = 0;
	int kkkk_debug = 0;

	void ToolpathGenerator::GenerateToolpathBasedonImageSpace()
	{
		BuildeImageSpace();

		int save_input_int_2 = input_int_2;
		input_int_2 = 100000000;

		clock_t start_time = clock();

		while (true)
		{
			image_space.Reset();
			image_space.CheckBoundary();
			image_space.UpdateDistanceFromBoundary();

			goon = false;
			over = false;
			spiral_affect_index = 0;
			iii_debug = 0;
			kkkk_debug = 0;

			int start_x, start_y;
			if (FindStartPixel(start_x, start_y))
			{
				SearchUnfilledRegionBoundary(image_space.pixels[start_x][start_y].related_boundary_x, image_space.pixels[start_x][start_y].related_boundary_y);
				InsertPixeltoSpiral(PixelIndex(start_x, start_y));

				if (save_input_int_2>=0)
				{
					input_int_2 = save_input_int_2;
				}

				SearchOnePath(true);

				if (!over)
				{
					break;
				}

				spirals.push_back(spiral);
				std::vector<PixelIndex>().swap(spiral);
				std::vector<PixelIndex>().swap(unfilled_region_boundary);
			}
			else
			{
				break;
			}
		}

		clock_t finish_time = clock();
		double m_time = (double)(finish_time - start_time) / CLOCKS_PER_SEC;
		std::cout << "Time: " << m_time << std::endl;
	}

	double ToolpathGenerator::SpiralPathDistance(int i, int j)
	{
		assert(i >= 0 && i<spiral.size() && j >= 0 && j<spiral.size());
		return image_space.PixelsDistance(spiral[i].x, spiral[i].y, spiral[j].x, spiral[j].y);
	}

	void ToolpathGenerator::SearchUnfilledRegionBoundary(int start_x, int start_y)
	{
		std::vector<PixelIndex> queue;

		int index_x = start_x;
		int index_y = start_y;

		queue.push_back(PixelIndex(index_x, index_y));

		while (queue.size() > 0)
		{
			index_x = queue[0].x;
			index_y = queue[0].y;
			queue.erase(queue.begin());
			unfilled_region_boundary.push_back(PixelIndex(index_x, index_y));

			image_space.pixels[index_x][index_y].boundary_index = unfilled_region_boundary.size() - 1;

			if (unfilled_region_boundary.size() == 1)
			{
				std::vector<PixelIndex> temp;
				for (int i = -1; i <= 1; i++)
				{
					for (int j = -1; j <= 1; j++)
					{
						if (!(i == 0 && j == 0) && image_space.check(index_x + i, index_y + j) && image_space.pixels[index_x + i][index_y + j].boundary)
						{
							temp.push_back(PixelIndex(index_x + i, index_y + j));
						}
					}
				}

				assert(temp.size()>0);

				//find the possible one
				int find_future_l = -1;
				bool find_b = false;

				for (int l = 0; l < temp.size() && !find_b; l++)
				{
					for (int i = -1; i <= 1 && !find_b; i++)
					{
						for (int j = -1; j <= 1 && !find_b; j++)
						{
							if (!(i == 0 && j == 0) && image_space.check(temp[l].x + i, temp[l].y + j) && image_space.pixels[temp[l].x + i][temp[l].y + j].boundary)
							{
								bool b = true;
								for (int k = 0; k < unfilled_region_boundary.size(); k++)
								{
									if (unfilled_region_boundary[k].x == (temp[l].x + i) && unfilled_region_boundary[k].y == (temp[l].y + j))
									{
										b = false;
										break;
									}
								}

								for (int k = 0; k < temp.size(); k++)
								{
									if (temp[k].x == (temp[l].x + i) && temp[k].y == (temp[l].y + j))
									{
										b = false;
										break;
									}
								}

								if (b)
								{
									find_future_l = l;
									find_b = true;
									break;
								}
							}
						}
					}
				}
				assert(find_b);

				for (int l = 0; l<temp.size(); l++)
				{
					if (l == find_future_l)
					{
						queue.push_back(temp[l]);
					}
					else
					{
						unfilled_region_boundary.push_back(temp[l]);
						image_space.pixels[temp[l].x][temp[l].y].boundary_index = unfilled_region_boundary.size() - 1;
					}
				}

				std::vector<PixelIndex>().swap(temp);
			}
			else
			{
				for (int i = -1; i <= 1; i++)
				{
					for (int j = -1; j <= 1; j++)
					{
						if (!(i == 0 && j == 0) && image_space.check(index_x + i, index_y + j) && image_space.pixels[index_x + i][index_y + j].boundary)
						{
							bool b = true;
							for (int k = 0; k < unfilled_region_boundary.size(); k++)
							{
								if (unfilled_region_boundary[k].x == (index_x + i) && unfilled_region_boundary[k].y == (index_y + j))
								{
									b = false;
									break;
								}
							}

							for (int k = 0; k < queue.size(); k++)
							{
								if (queue[k].x == (index_x + i) && queue[k].y == (index_y + j))
								{
									b = false;
									break;
								}
							}

							if (b)
							{
								queue.push_back(PixelIndex(index_x + i, index_y + j));
							}
						}
					}
				}
			}
		}

		std::vector<PixelIndex>().swap(queue);

		//post progress
		int t = 0;

		for (int i = 1; i < unfilled_region_boundary.size()-1; i++)
		{
			if (!(image_space.PixelsNeighbour(unfilled_region_boundary[i].x, unfilled_region_boundary[i].y, unfilled_region_boundary[i - 1].x, unfilled_region_boundary[i - 1].y) &&
				image_space.PixelsNeighbour(unfilled_region_boundary[i].x, unfilled_region_boundary[i].y, unfilled_region_boundary[i + 1].x, unfilled_region_boundary[i + 1].y) &&
				!image_space.PixelsNeighbour(unfilled_region_boundary[i - 1].x, unfilled_region_boundary[i - 1].y, unfilled_region_boundary[i + 1].x, unfilled_region_boundary[i + 1].y)))
			{
				image_space.pixels[unfilled_region_boundary[i].x][unfilled_region_boundary[i].y].boundary_index = t;
			}
			else
			{
				image_space.pixels[unfilled_region_boundary[i].x][unfilled_region_boundary[i].y].boundary_index = i;
				t=i;
			}
		}
	}

	void ToolpathGenerator::InsertPixeltoSpiral(PixelIndex pi)
	{
		
		if (spiral.size()>0 && image_space.pixels[pi.x][pi.y].related_boundary_x >= 0 && image_space.pixels[pi.x][pi.y].related_boundary_y >= 0)
		{
			if (image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_x >= 0 && image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_y>=0)
			{
				int int1 = image_space.pixels[image_space.pixels[pi.x][pi.y].related_boundary_x][image_space.pixels[pi.x][pi.y].related_boundary_y].boundary_index;
				int int2 = image_space.pixels[image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_x][image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_y].boundary_index;
			
				if (int2 - int1 > 0)
				{
					//std::cerr << "Boundary index: " << int2 << " " << int1 << " :----------" << std::endl;
					spiral_direction = false;
				}

				if (int2 == int1)
				{
					//std::cerr << "Boundary index: " << int2 << " " << int1 << " :=========" << std::endl;
				}

				if (int2 - int1 < 0)
				{
					//std::cerr << "Boundary index: " << int2 << " " << int1 << " :+++++++++" << std::endl;
					spiral_direction = true;
				}
			}
		}

		image_space.pixels[pi.x][pi.y].centeraxis = true;
		spiral.push_back(pi);
		
		for (int i = pi.x - 1.0*toolpath_size / image_space.pixel_size; i <= pi.x + 1.0*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = pi.y - 1.0*toolpath_size / image_space.pixel_size; j <= pi.y + 1.0*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j) && !image_space.pixels[i][j].filled)
				{
					double doubletemp = image_space.PixelsDistance(i, j, pi.x, pi.y);
					if (image_space.PixelsDistance(i, j, pi.x, pi.y)<=toolpath_size / 2.0+0.00001)
					{
						image_space.pixels[i][j].filled = true;
						image_space.pixels[i][j].filled_spiral_index = spiral.size() - 1;
					}
				}
			}
		}

		while (SpiralPathDistance(spiral.size() - 1, spiral_affect_index) > (2.0*toolpath_size + image_space.pixel_size))
		{
			for (int i = spiral[spiral_affect_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral[spiral_affect_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
			{
				for (int j = spiral[spiral_affect_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral[spiral_affect_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
				{
					if (image_space.check(i, j))
					{
						image_space.pixels[i][j].save_distance = image_space.pixels[i][j].distance;

						image_space.pixels[i][j].save_related_boundary_index = image_space.pixels[i][j].related_boundary_index;
						image_space.pixels[i][j].save_related_boundary_x = image_space.pixels[i][j].related_boundary_x;
						image_space.pixels[i][j].save_related_boundary_y = image_space.pixels[i][j].related_boundary_y;
					}
				}
			}

			for (int i = spiral[spiral_affect_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral[spiral_affect_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
			{
				for (int j = spiral[spiral_affect_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral[spiral_affect_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
				{
					if (image_space.check(i, j))
					{
						if (image_space.pixels[i][j].distance > image_space.PixelsDistance(i, j, spiral[spiral_affect_index].x, spiral[spiral_affect_index].y)-toolpath_size/2.0)
						{
							image_space.pixels[i][j].distance = image_space.PixelsDistance(i, j, spiral[spiral_affect_index].x, spiral[spiral_affect_index].y) - toolpath_size / 2.0;
							image_space.pixels[i][j].related_boundary_index = spiral_affect_index;
							image_space.pixels[i][j].related_boundary_x = spiral[spiral_affect_index].x;
							image_space.pixels[i][j].related_boundary_y = spiral[spiral_affect_index].y;
						}
					}
				}
			}

			for (int i = spiral[spiral_affect_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral[spiral_affect_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
			{
				for (int j = spiral[spiral_affect_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral[spiral_affect_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
				{
					if (image_space.check(i, j))
					{
						if (image_space.pixels[i][j].distance < 1.0*toolpath_size + image_space.pixel_size&&
							!(image_space.pixels[i][j].related_boundary_x == image_space.pixels[i][j].save_related_boundary_x&&
							image_space.pixels[i][j].related_boundary_y == image_space.pixels[i][j].save_related_boundary_y))
						{
							image_space.pixels[i][j].bool_t = true;
						}
						else
						{
							image_space.pixels[i][j].distance = image_space.pixels[i][j].save_distance;

							image_space.pixels[i][j].related_boundary_index = image_space.pixels[i][j].save_related_boundary_index;
							image_space.pixels[i][j].related_boundary_x = image_space.pixels[i][j].save_related_boundary_x;
							image_space.pixels[i][j].related_boundary_y = image_space.pixels[i][j].save_related_boundary_y;
						}
					}
				}
			}
			spiral_affect_index++;
		}
	}

	void ToolpathGenerator::SearchOnePath(bool b)
	{
		int jjj_debug = 0;
		
		do
		{
			if (over)
				break;

			jjj_debug++;
			kkkk_debug++;
			iii_debug++;

			std::cout << kkkk_debug << std::endl;

			goon = false;
			double min_d = MAXDOUBLE;
			int index_x = -1;
			int index_y = -1;

			std::vector<PixelIndex> perfect_index;
			std::vector<PixelIndex> all_index;
			std::vector<double> all_index_error;
			std::vector<double> all_index_distance;


			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					if (image_space.check(spiral[spiral.size() - 1].x + i, spiral[spiral.size() - 1].y + j))
					{
						Pixel &p = image_space.pixels[spiral[spiral.size() - 1].x + i][spiral[spiral.size() - 1].y + j];

						int current_index_x = spiral[spiral.size() - 1].x + i;
						int current_index_y = spiral[spiral.size() - 1].y + j;

						if (p.inside&&!p.centeraxis)
						{
							bool use_this_pixel = true;

							for (int x = current_index_x - toolpath_size / image_space.pixel_size; x <= current_index_x + toolpath_size / image_space.pixel_size; x++)
							{
								for (int y = current_index_y - toolpath_size / image_space.pixel_size; y <= current_index_y + 1.0*toolpath_size / image_space.pixel_size; y++)
								{
									if (image_space.check(x, y) && image_space.PixelsDistance(current_index_x, current_index_y, x, y) < toolpath_size / 2.0 - image_space.pixel_size&&use_this_pixel)
									{
										if (image_space.pixels[x][y].filled&&image_space.pixels[x][y].filled_spiral_index>=0&&SpiralPathDistance(spiral.size() - 1, image_space.pixels[x][y].filled_spiral_index)>toolpath_size + image_space.pixel_size)
										{
											use_this_pixel = false;
											break;
										}
									}
								}
							}

							if (spiral_affect_index > 0)
							{
								if (spiral.size()>0 && image_space.pixels[current_index_x][current_index_y].related_boundary_x >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_y >= 0)
								{
									if (image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_x >= 0 && image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_y >= 0)
									{
										int int1 = image_space.pixels[image_space.pixels[current_index_x][current_index_y].related_boundary_x][image_space.pixels[current_index_x][current_index_y].related_boundary_y].boundary_index;
										int int2 = image_space.pixels[image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_x][image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_y].boundary_index;

										if (int1 >= 0 && int2 >= 0)
										{
											if (int2 - int1 > 0)
											{
												if (spiral_direction)
												{
													use_this_pixel = false;
												}
											}

											if (int2 - int1 < 0)
											{
												if (!spiral_direction)
												{
													use_this_pixel = false;
												}
											}
										}
									}
								}
							}

							if (spiral.size()>0 && image_space.pixels[current_index_x][current_index_y].related_boundary_x >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_y >= 0)
							{
								if (image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_x >= 0 && image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_y >= 0)
								{
									if (image_space.pixels[current_index_x][current_index_y].related_boundary_index < image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_index)
									{
										use_this_pixel = false;
									}
								}
							}

							if (use_this_pixel)
							{
								if (abs(p.distance - toolpath_size / 2.0) < 0.0001)
								{
									perfect_index.push_back(PixelIndex(spiral[spiral.size() - 1].x + i, spiral[spiral.size() - 1].y + j));
								}

								all_index.push_back(PixelIndex(i, j));
								all_index_error.push_back(abs(p.distance - toolpath_size / 2.0));
								all_index_distance.push_back(p.distance);

								if (abs(p.distance - toolpath_size / 2.0) < min_d)
								{
									min_d = abs(p.distance - toolpath_size / 2.0);
									index_x = spiral[spiral.size() - 1].x + i;
									index_y = spiral[spiral.size() - 1].y + j;
									goon = true;
								}
							}
						}
					}
				}
			}

			if (perfect_index.size() > 1)
			{
				int spiral_index_min = (int)MAXDOUBLE;
				int spiral_index;
				for (int i = 0; i < perfect_index.size(); i++)
				{
					if (image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index>=image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_index)
					{
						if (spiral_index_min > image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index)
						{
							spiral_index_min = image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index;
							spiral_index = i;
						}
					}
				}

				index_x = perfect_index[spiral_index].x;
				index_y = perfect_index[spiral_index].y;
				goon = true;
			}
			else
			{
				int index_x_x = index_x;
				int index_y_y = index_y;

				if (spiral.size() >= 2)
				{

					double max_distance_temp = image_space.PixelsDistance(index_x, index_y, spiral[spiral.size() - 2].x, spiral[spiral.size() - 2].y);

					for (int i = 0; i < all_index_error.size(); i++)
					{
						if (!(all_index[i].x == (index_x - spiral[spiral.size() - 1].x) && all_index[i].y == (index_y - spiral[spiral.size() - 1].y)))
						{
							if (abs(all_index_error[i] - min_d) < 0.0001)
							{
								double distance_temp = image_space.PixelsDistance(spiral[spiral.size() - 1].x + all_index[i].x, spiral[spiral.size() - 1].y + all_index[i].y, spiral[spiral.size() - 2].x, spiral[spiral.size() - 2].y);

								if (image_space.pixels[index_x_x][index_y_y].related_boundary_index == image_space.pixels[spiral[spiral.size() - 1].x + all_index[i].x][spiral[spiral.size() - 1].y + all_index[i].y].related_boundary_index)
								{
									if (max_distance_temp < distance_temp)
									{
										max_distance_temp = distance_temp;

										index_x_x = spiral[spiral.size() - 1].x + all_index[i].x;
										index_y_y = spiral[spiral.size() - 1].y + all_index[i].y;
									}
								}
								else
								{
									if (image_space.pixels[index_x_x][index_y_y].related_boundary_index<image_space.pixels[spiral[spiral.size() - 1].x + all_index[i].x][spiral[spiral.size() - 1].y + all_index[i].y].related_boundary_index)
									{
										max_distance_temp = distance_temp;
										index_x_x = spiral[spiral.size() - 1].x + all_index[i].x;
										index_y_y = spiral[spiral.size() - 1].y + all_index[i].y;
									}
								}
							}
						}
					}
				}

				index_x = index_x_x;
				index_y = index_y_y;
			}
			
			std::vector<PixelIndex>().swap(perfect_index);
			std::vector<PixelIndex>().swap(all_index);
			std::vector<double>().swap(all_index_error);
			std::vector<double>().swap(all_index_distance);

			if (goon)
			{
				goon = false;
				for (int x = index_x - toolpath_size / image_space.pixel_size; x <= index_x + toolpath_size / image_space.pixel_size; x++)
				{
					for (int y = index_y - toolpath_size / image_space.pixel_size; y <= index_y +toolpath_size / image_space.pixel_size; y++)
					{
						if (image_space.check(x, y) && image_space.PixelsDistance(index_x, index_y, x, y) <= toolpath_size / 2.0 + 0.00001 && !goon)
						{
							if (!image_space.pixels[x][y].filled)
							{
								goon = true;
								break;
							}
						}
					}
				}
			}

			if (goon)
			{
				InsertPixeltoSpiral(PixelIndex(index_x, index_y));
			}
			else
			{
				over = true;
				std::cout << "Over..." << std::endl;
				break;
			}

			if (!b)
			{
				iii_debug=0;
			}

			if (!b&&jjj_debug==input_int_1)
			{
				break;
			}

		//} while (true);
		} while (iii_debug<input_int_2);

	}

	bool ToolpathGenerator::FindStartPixel(int &index_x, int &index_y)
	{
		double start_pixel_d = MAXDOUBLE;

		for (int i = 0; i < image_space.pixels.size(); i++)
		{
			for (int j = 0; j < image_space.pixels[i].size(); j++)
			{
				if (image_space.pixels[i][j].inside&& !image_space.pixels[i][j].filled)
				{
					if (abs(image_space.pixels[i][j].distance - toolpath_size / 2.0) < start_pixel_d)
					{
						start_pixel_d = abs(image_space.pixels[i][j].distance - toolpath_size / 2.0);
						index_x = i;
						index_y = j;
					}
				}
			}
		}

		for (int i = 0; i < image_space.pixels.size() && !goon; i++)
		{
			for (int j = 0; j < image_space.pixels[i].size() && !goon; j++)
			{
				if (image_space.pixels[i][j].inside&& !image_space.pixels[i][j].filled)
				{
					if (image_space.pixels[i][j].distance >= toolpath_size / 2.0)
					{
						return true;
					}
				}
			}
		}

		return false;
	}

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

}