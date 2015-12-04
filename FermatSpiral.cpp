#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"


namespace hpcg {

	void ToolpathGenerator::FermatSpiral()
	{
		GenerateToolpath();
		BuildeImageSpace();
		ChooseEntryExitPoints();
	}

	int entry_number = 0;
	int exit_number = 0;

	void ToolpathGenerator::GetOnePointFromOffset(double d, int &index_x, int &index_y)
	{
		assert(d>=0.0&&d<=1.0);

		double length = 0.0;
		for (int i = 0; i < toolpath[0].size(); i++)
		{
			length += sqrt((double)CGAL::squared_distance(Point_2(toolpath[0][i].x, toolpath[0][i].y), Point_2(toolpath[0][(i + 1) % toolpath[0].size()].x, toolpath[0][(i + 1) % toolpath[0].size()].y)));
		}

		double total_length = length;
		length = 0.0;

		for (int i = 0; i < toolpath[0].size(); i++)
		{
			double l = sqrt((double)CGAL::squared_distance(Point_2(toolpath[0][i].x, toolpath[0][i].y), Point_2(toolpath[0][(i + 1) % toolpath[0].size()].x, toolpath[0][(i + 1) % toolpath[0].size()].y)));

			if (d*total_length >= length&&d*total_length <= length + l)
			{
				double ll = (d - length / total_length)*total_length/l;
				double x = toolpath[0][i].x + (toolpath[0][(i + 1) % toolpath[0].size()].x - toolpath[0][i].x)*ll;
				double y = toolpath[0][i].y + (toolpath[0][(i + 1) % toolpath[0].size()].y - toolpath[0][i].y)*ll;

				image_space.FindClosestPoint(x, y, index_x, index_y);
				break;
			}
			length += l;
		}

		double min_d = abs(image_space.pixels[index_x][index_y].distance-toolpath_size/2.0);

		int i_int = 0;
		int j_int = 0;
		for (int i = -2; i <= 2; i++)
		{
			for (int j = -2; j <= 2; j++)
			{
				if (!(i == 0 && j == 0) && image_space.check(index_x + i, index_y + j))
				{
					if (abs(image_space.pixels[index_x + i][index_y + j].distance - toolpath_size / 2.0) < min_d)
					{
						min_d = abs(image_space.pixels[index_x + i][index_y + j].distance - toolpath_size / 2.0);
						i_int = i;
						j_int = j;
					}
				}
			}
		}

		index_x = index_x + i_int;
		index_y = index_y + j_int;
	}
	
	void ToolpathGenerator::ConnectTwoPixels(PixelIndex p1, PixelIndex p2, std::vector<PixelIndex> &insertpixels)
	{
		int x = p1.x;
		int y = p1.y;

		Line_2 line(Point_2(image_space.pixels[p2.x][p2.y].center[0], image_space.pixels[p2.x][p2.y].center[1]), Point_2(image_space.pixels[p1.x][p1.y].center[0], image_space.pixels[p1.x][p1.y].center[1]));

		double max_d = image_space.PixelsDistance(p1.x, p1.y, p2.x, p2.y);

		insertpixels.push_back(p1);
		do
		{
			double min_d = MAXDOUBLE;
			int index_x = -1;
			int index_y = -1;
			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					if (image_space.check(x + i, y + j) && !(i == 0 && j == 0) && !((x + i) == p1.x && (y + j) == p1.y))
					{
						if (image_space.PixelsDistance(x + i, y + j, p2.x, p2.y) <= max_d)
						{
							Pixel &p = image_space.pixels[x + i][y + j];

							if (insertpixels.size() > 1)
							{
								if ((x + i) == insertpixels[insertpixels.size() - 2].x && (y + j) == insertpixels[insertpixels.size() - 2].y)
								{
									continue;
								}
							}
							double d = sqrt(CGAL::squared_distance(Point_2(p.center[0], p.center[1]), line));
							if (min_d > d)
							{
								min_d = d;
								index_x = x + i;
								index_y = y + j;
							}
						}
					}
				}
			}

			if (index_x >= 0 && index_y >= 0)
			{
				if (index_x == p2.x&&index_y == p2.y)
				{
					break;
				}
				insertpixels.push_back(PixelIndex(index_x, index_y));
			}
			x = index_x;
			y = index_y;

		} while (true);

		insertpixels.erase(insertpixels.begin());
	}

	void ToolpathGenerator::InsertPixeltoSpiralEntry(PixelIndex pi)
	{
		if (!first_two_step_bool&&spiral_entry.size()>0 && image_space.pixels[pi.x][pi.y].related_boundary_x >= 0 && image_space.pixels[pi.x][pi.y].related_boundary_y >= 0)
		{
			if (image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x >= 0 && image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y >= 0)
			{
				int int1 = image_space.pixels[image_space.pixels[pi.x][pi.y].related_boundary_x][image_space.pixels[pi.x][pi.y].related_boundary_y].boundary_index;
				int int2 = image_space.pixels[image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x][image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y].boundary_index;

				if (int2 - int1 > 0)
				{
					spiral_direction = false;
				}

				if (int2 - int1 < 0)
				{
					spiral_direction = true;
				}
			}
		}

		image_space.pixels[pi.x][pi.y].centeraxis = true;
		spiral_entry.push_back(pi);

		for (int i = pi.x - 1.0*toolpath_size / image_space.pixel_size; i <= pi.x + 1.0*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = pi.y - 1.0*toolpath_size / image_space.pixel_size; j <= pi.y + 1.0*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j) && !image_space.pixels[i][j].filled)
				{
					if (image_space.PixelsDistance(i, j, pi.x, pi.y) <= toolpath_size / 2.0 + 0.00001)
					{
						image_space.pixels[i][j].filled = true;
						image_space.pixels[i][j].filled_entry_exit = true;
						image_space.pixels[i][j].filled_spiral_index = spiral_entry.size() - 1;
					}
				}
			}
		}



		int current_index = spiral_entry.size() - 1;

		for (spiral_entry_affect_index = current_index; spiral_entry_affect_index >= 0; spiral_entry_affect_index--)
		{
			if (image_space.PixelsDistance(spiral_entry[spiral_entry.size() - 1].x, spiral_entry[spiral_entry.size() - 1].y, spiral_entry[spiral_entry_affect_index].x, spiral_entry[spiral_entry_affect_index].y) > (2.0*toolpath_size + image_space.pixel_size))
			{
				break;
			}
		}

		for (int i = spiral_entry[current_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral_entry[current_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = spiral_entry[current_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral_entry[current_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j))
				{
					image_space.pixels[i][j].save_distance = image_space.pixels[i][j].distance_entry;
					image_space.pixels[i][j].save_related_boundary_index = image_space.pixels[i][j].related_boundary_index_entry;
					image_space.pixels[i][j].save_related_boundary_x = image_space.pixels[i][j].related_boundary_x_entry;
					image_space.pixels[i][j].save_related_boundary_y = image_space.pixels[i][j].related_boundary_y_entry;
				}
			}
		}

		for (int i = spiral_entry[current_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral_entry[current_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = spiral_entry[current_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral_entry[current_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j))
				{
					double d = image_space.PixelsDistance(i, j, spiral_entry[current_index].x, spiral_entry[current_index].y) - toolpath_size / 2.0;

					if (image_space.pixels[i][j].distance_entry > d&&abs(image_space.pixels[i][j].distance_entry - d)>0.00001)
					{
						image_space.pixels[i][j].distance_entry = d;
						image_space.pixels[i][j].related_boundary_index_entry = current_index;
						image_space.pixels[i][j].related_boundary_x_entry = spiral_entry[current_index].x;
						image_space.pixels[i][j].related_boundary_y_entry = spiral_entry[current_index].y;
					}
				}
			}
		}

		for (int i = spiral_entry[current_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral_entry[current_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = spiral_entry[current_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral_entry[current_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j))
				{
					if (image_space.pixels[i][j].distance_entry < 1.0*toolpath_size + image_space.pixel_size&&
						!(image_space.pixels[i][j].related_boundary_x_entry == image_space.pixels[i][j].save_related_boundary_x&&
						image_space.pixels[i][j].related_boundary_y_entry == image_space.pixels[i][j].save_related_boundary_y))
					{
						image_space.pixels[i][j].bool_t = true;
						image_space.pixels[i][j].bool_t_entry_exit = true;
						image_space.pixels[i][j].bool_t_entry = true;
					}
					else
					{
						image_space.pixels[i][j].distance_entry = image_space.pixels[i][j].save_distance;
						image_space.pixels[i][j].related_boundary_index_entry = image_space.pixels[i][j].save_related_boundary_index;
						image_space.pixels[i][j].related_boundary_x_entry = image_space.pixels[i][j].save_related_boundary_x;
						image_space.pixels[i][j].related_boundary_y_entry = image_space.pixels[i][j].save_related_boundary_y;
					}
				}
			}
		}
	}

	void ToolpathGenerator::InsertPixeltoSpiralExit(PixelIndex pi)
	{
		image_space.pixels[pi.x][pi.y].centeraxis = true;
		spiral_exit.push_back(pi);

		for (int i = pi.x - 1.0*toolpath_size / image_space.pixel_size; i <= pi.x + 1.0*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = pi.y - 1.0*toolpath_size / image_space.pixel_size; j <= pi.y + 1.0*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j) && !image_space.pixels[i][j].filled)
				{
					if (image_space.PixelsDistance(i, j, pi.x, pi.y) <= toolpath_size / 2.0+0.00001)
					{
						image_space.pixels[i][j].filled = true;
						image_space.pixels[i][j].filled_entry_exit = false;
						image_space.pixels[i][j].filled_spiral_index = spiral_exit.size() - 1;
					}
				}
			}
		}

		int current_index = spiral_exit.size() - 1;

		for (spiral_exit_affect_index = current_index; spiral_exit_affect_index >= 0; spiral_exit_affect_index--)
		{
			if (image_space.PixelsDistance(spiral_exit[spiral_exit.size() - 1].x, spiral_exit[spiral_exit.size() - 1].y, spiral_exit[spiral_exit_affect_index].x, spiral_exit[spiral_exit_affect_index].y) > (2.0*toolpath_size + image_space.pixel_size))
			{
				break;
			}
		}

		for (int i = spiral_exit[current_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral_exit[current_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = spiral_exit[current_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral_exit[current_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j))
				{
					image_space.pixels[i][j].save_distance = image_space.pixels[i][j].distance_exit;
					image_space.pixels[i][j].save_related_boundary_index = image_space.pixels[i][j].related_boundary_index_exit;
					image_space.pixels[i][j].save_related_boundary_x = image_space.pixels[i][j].related_boundary_x_exit;
					image_space.pixels[i][j].save_related_boundary_y = image_space.pixels[i][j].related_boundary_y_exit;
				}
			}
		}

		for (int i = spiral_exit[current_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral_exit[current_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = spiral_exit[current_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral_exit[current_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j))
				{
					double d = image_space.PixelsDistance(i, j, spiral_exit[current_index].x, spiral_exit[current_index].y) - toolpath_size / 2.0;

					if (image_space.pixels[i][j].distance_exit > d&&abs(image_space.pixels[i][j].distance_exit - d)>0.00001)
					{
						image_space.pixels[i][j].distance_exit = d;
						image_space.pixels[i][j].related_boundary_index_exit = current_index;
						image_space.pixels[i][j].related_boundary_x_exit = spiral_exit[current_index].x;
						image_space.pixels[i][j].related_boundary_y_exit = spiral_exit[current_index].y;
					}
				}
			}
		}

		for (int i = spiral_exit[current_index].x - 1.5*toolpath_size / image_space.pixel_size; i <= spiral_exit[current_index].x + 1.5*toolpath_size / image_space.pixel_size; i++)
		{
			for (int j = spiral_exit[current_index].y - 1.5*toolpath_size / image_space.pixel_size; j <= spiral_exit[current_index].y + 1.5*toolpath_size / image_space.pixel_size; j++)
			{
				if (image_space.check(i, j))
				{
					if (image_space.pixels[i][j].distance_exit < 1.0*toolpath_size + image_space.pixel_size&&
						!(image_space.pixels[i][j].related_boundary_x_exit == image_space.pixels[i][j].save_related_boundary_x&&
						image_space.pixels[i][j].related_boundary_y_exit == image_space.pixels[i][j].save_related_boundary_y))
					{
						image_space.pixels[i][j].bool_t = true;
						image_space.pixels[i][j].bool_t_entry_exit = false;
						image_space.pixels[i][j].bool_t_exit = true;
					}
					else
					{
						image_space.pixels[i][j].distance_exit = image_space.pixels[i][j].save_distance;
						image_space.pixels[i][j].related_boundary_index_exit = image_space.pixels[i][j].save_related_boundary_index;
						image_space.pixels[i][j].related_boundary_x_exit = image_space.pixels[i][j].save_related_boundary_x;
						image_space.pixels[i][j].related_boundary_y_exit = image_space.pixels[i][j].save_related_boundary_y;
					}
				}
			}
		}
	}

	int iter_exit = 0;
	int iter_entry = 0;

	void ToolpathGenerator::SearchOneStepforSpiralExit(bool b)
	{
		iter_entry = 0;
		bool goon = true;

		do
		{
			std::cerr << iter_exit << std::endl;
			iter_exit++;

			double min_d = MAXDOUBLE;
			int index_x = -1;
			int index_y = -1;

			std::vector<PixelIndex> perfect_index;
			std::vector<PixelIndex> all_index;
			std::vector<double> all_index_error;
			std::vector<double> all_index_distance;
			std::vector<Pixel> all_pixels;
		
			int c_x = spiral_exit[spiral_exit.size() - 1].x;
			int c_y = spiral_exit[spiral_exit.size() - 1].y;

			bool goon_check = true;
			goon = false;

			Polygon_2 one_contour;
			bool one_contour_bool = false;
			if (spiral_exit.size() > 1)
			{

				int related_boundary_x1 = image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_x_entry;
				int related_boundary_y1 = image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_y_entry;
				
				int related_boundary_x2 = image_space.pixels[spiral_exit[spiral_exit.size() - 2].x][spiral_exit[spiral_exit.size() - 2].y].related_boundary_x_entry;
				int related_boundary_y2 = image_space.pixels[spiral_exit[spiral_exit.size() - 2].x][spiral_exit[spiral_exit.size() - 2].y].related_boundary_y_entry;

				if (related_boundary_x1 >= 0 && related_boundary_y1 >= 0 && related_boundary_x2 >= 0 && related_boundary_y2 >= 0)
				{
					one_contour.push_back(Point_2(image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].center.x, image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].center.y));
					one_contour.push_back(Point_2(image_space.pixels[spiral_exit[spiral_exit.size() - 2].x][spiral_exit[spiral_exit.size() - 2].y].center.x, image_space.pixels[spiral_exit[spiral_exit.size() - 2].x][spiral_exit[spiral_exit.size() - 2].y].center.y));

					if (related_boundary_x1 == related_boundary_x2&&related_boundary_y1 == related_boundary_y2)
					{
						one_contour.push_back(Point_2(image_space.pixels[related_boundary_x2][related_boundary_y2].center.x, image_space.pixels[related_boundary_x2][related_boundary_y2].center.y));
					}
					else
					{
						one_contour.push_back(Point_2(image_space.pixels[related_boundary_x2][related_boundary_y2].center.x, image_space.pixels[related_boundary_x2][related_boundary_y2].center.y));
						one_contour.push_back(Point_2(image_space.pixels[related_boundary_x1][related_boundary_y1].center.x, image_space.pixels[related_boundary_x1][related_boundary_y1].center.y));
					}
					one_contour_bool = true;
				}
			}

			for (int i = -1; i <= 1 && goon_check; i++)
			{
				for (int j = -1; j <= 1 && goon_check; j++)
				{
					if (image_space.check(spiral_exit[spiral_exit.size() - 1].x + i, spiral_exit[spiral_exit.size() - 1].y + j))
					{
						Pixel &p = image_space.pixels[spiral_exit[spiral_exit.size() - 1].x + i][spiral_exit[spiral_exit.size() - 1].y + j];

						int current_index_x = spiral_exit[spiral_exit.size() - 1].x + i;
						int current_index_y = spiral_exit[spiral_exit.size() - 1].y + j;

						if (p.inside&&!p.centeraxis)
						{
							bool use_this_pixel = true;

							if (image_space.pixels[current_index_x][current_index_y].failure_nb >= 2)
							{
								use_this_pixel = false;
							}

							if (spiral_exit.size() == 1||first_two_step_bool)
							{
								for (int x = current_index_x - toolpath_size / image_space.pixel_size; x <= current_index_x + toolpath_size / image_space.pixel_size; x++)
								{
									for (int y = current_index_y - toolpath_size / image_space.pixel_size; y <= current_index_y + 1.0*toolpath_size / image_space.pixel_size; y++)
									{
										if (image_space.check(x, y) && image_space.PixelsDistance(current_index_x, current_index_y, x, y) <= toolpath_size / 2.0 + 0.00001&& use_this_pixel)
										{
											if (image_space.pixels[x][y].filled&&image_space.pixels[x][y].filled_entry_exit)
											{
												use_this_pixel = false;
												break;
											}
										}
									}
								}
							}
				
							//if (false)
							if (spiral_exit.size() > 1)
							{
								int related_boundary_x_entry = image_space.pixels[current_index_x][current_index_y].related_boundary_x_entry;
								int related_boundary_y_entry = image_space.pixels[current_index_x][current_index_y].related_boundary_y_entry;

								if (image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_x_entry == related_boundary_x_entry&&
									image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_y_entry == related_boundary_y_entry&&
									image_space.pixels[spiral_exit[spiral_exit.size() - 2].x][spiral_exit[spiral_exit.size() - 2].y].related_boundary_x_entry == related_boundary_x_entry&&
									image_space.pixels[spiral_exit[spiral_exit.size() - 2].x][spiral_exit[spiral_exit.size() - 2].y].related_boundary_y_entry == related_boundary_y_entry)
								{
									Vector2d v0 = image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].center - image_space.pixels[related_boundary_x_entry][related_boundary_y_entry].center;
									Vector2d v1 = image_space.pixels[spiral_exit[spiral_exit.size() - 2].x][spiral_exit[spiral_exit.size() - 2].y].center - image_space.pixels[related_boundary_x_entry][related_boundary_y_entry].center;
									Vector2d v2 = image_space.pixels[current_index_x][current_index_y].center - image_space.pixels[related_boundary_x_entry][related_boundary_y_entry].center;

									if ((v1[0] * v0[1] - v1[1] * v0[0])*(v2[0] * v0[1] - v2[1] * v0[0]) > 0)
									{
										use_this_pixel = false;
									}
								}
							}
						
							if (one_contour_bool&&one_contour.is_simple())
							{
								if (one_contour.bounded_side(Point_2(image_space.pixels[current_index_x][current_index_y].center.x, image_space.pixels[current_index_x][current_index_y].center.y)) == CGAL::ON_UNBOUNDED_SIDE)
								{
									int related_boundary_x_entry = image_space.pixels[current_index_x][current_index_y].related_boundary_x_entry;
									int related_boundary_y_entry = image_space.pixels[current_index_x][current_index_y].related_boundary_y_entry;

									int related_boundary_x1 = image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_x_entry;
									int related_boundary_y1 = image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_y_entry;

									if (related_boundary_x_entry == related_boundary_x1&&related_boundary_y_entry == related_boundary_y1)
									{
										Point_2 p0(image_space.pixels[current_index_x][current_index_y].center.x, image_space.pixels[current_index_x][current_index_y].center.y);
										Point_2 p1(image_space.pixels[related_boundary_x_entry][related_boundary_y_entry].center.x, image_space.pixels[related_boundary_x_entry][related_boundary_y_entry].center.y);

										Segment_2 s0(p0, p1);

										bool center_bool = false;

										if (one_contour.bounded_side(Point_2(image_space.pixels[current_index_x][current_index_y].center.x, image_space.pixels[current_index_x][current_index_y].center.y)) != CGAL::ON_BOUNDARY)
										{
											center_bool = true;
										}

										bool related_bool = false;

										std::vector<Point_2> ignore_ps;

										int ver_nb = 0;
										for (Polygon_2::Vertex_iterator ver_iter = one_contour.vertices_begin(); ver_iter != one_contour.vertices_end(); ver_iter++)
										{
											ver_nb++;
										}
										for (Polygon_2::Vertex_iterator ver_iter = one_contour.vertices_begin(); ver_iter != one_contour.vertices_end(); ver_iter++)
										{

											double  xxx1 = ver_iter->x();
											double  xxx2 = ver_iter->y();

											Polygon_2::Vertex_iterator pre_ver_iter0 = one_contour.vertices_begin();
											Polygon_2::Vertex_iterator pre_ver_iter1 = one_contour.vertices_begin() + 1;
											Polygon_2::Vertex_iterator pre_ver_iter2 = one_contour.vertices_begin() + 2;

											if (abs(ver_iter->x() - p1.x()) < 0.0000001&&abs(ver_iter->y() - p1.y()) < 0.0000001)
											{
												related_bool = true;
												ignore_ps.push_back(Point_2(p1.x(), p1.y()));

												if (ver_iter == one_contour.vertices_begin())
												{
													Polygon_2::Vertex_iterator pre_ver_iter = one_contour.vertices_begin() + ver_nb - 1;
													ignore_ps.push_back(Point_2(pre_ver_iter->x(), pre_ver_iter->y()));
												}
												else
												{
													Polygon_2::Vertex_iterator pre_ver_iter = ver_iter - 1;
													ignore_ps.push_back(Point_2(pre_ver_iter->x(), pre_ver_iter->y()));
												}

												Polygon_2::Vertex_iterator next_ver_iter = ver_iter + 1;
												if (next_ver_iter == one_contour.vertices_end())
												{
													ignore_ps.push_back(Point_2(one_contour.vertices_begin()->x(), one_contour.vertices_begin()->y()));
												}
												else
												{
													ignore_ps.push_back(Point_2(next_ver_iter->x(), next_ver_iter->y()));
												}
												break;
											}
										}

										for (Polygon_2::Edge_const_iterator edge_iter = one_contour.edges_begin(); edge_iter != one_contour.edges_end(); edge_iter++)
										{
											Segment_2 s1(edge_iter->source(), edge_iter->target());
											CGAL::Object result = CGAL::intersection(s0, s1);
											if (const Point_2 *ipoint = CGAL::object_cast<Point_2 >(&result)) {
												if (related_bool)
												{
													if (ignore_ps.size() == 3)
													{
														if (!(abs(ignore_ps[0].x() - ipoint->x()) < 0.0000001&&abs(ignore_ps[0].y() - ipoint->y()) < 0.0000001) &&
															!(abs(ignore_ps[1].x() - ipoint->x()) < 0.0000001&&abs(ignore_ps[1].y() - ipoint->y()) < 0.0000001) &&
															!(abs(ignore_ps[2].x() - ipoint->x()) < 0.0000001&&abs(ignore_ps[2].y() - ipoint->y()) < 0.0000001))
														{
															use_this_pixel = false;
															break;
														}
													}
												}
												else
												{
													use_this_pixel = false;
													break;
												}
											}
											else
												if (const Segment_2 *iseg = CGAL::object_cast<Segment_2>(&result)) {

												}
												else {

												}
										}
										std::vector<Point_2>().swap(ignore_ps);
									}
								}
								else
								{
									use_this_pixel = false;
								}
							}

							if (spiral_exit_affect_index > 0 && !first_two_step_bool)
							{
								if (spiral_exit.size()>0 && image_space.pixels[current_index_x][current_index_y].related_boundary_x_entry >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_y_entry >= 0)
								{
									if (image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_x_entry >= 0 && image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_y_entry >= 0)
									{
										int int1 = image_space.pixels[image_space.pixels[current_index_x][current_index_y].related_boundary_x_entry][image_space.pixels[current_index_x][current_index_y].related_boundary_y_entry].boundary_index;
										int int2 = image_space.pixels[image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_x_entry][image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_y_entry].boundary_index;

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

							if (first_two_step_bool)
							for (int x = current_index_x - toolpath_size / image_space.pixel_size; x <= current_index_x + toolpath_size / image_space.pixel_size; x++)
							{
								for (int y = current_index_y - toolpath_size / image_space.pixel_size; y <= current_index_y + 1.0*toolpath_size / image_space.pixel_size; y++)
								{
									if (image_space.check(x, y) && image_space.PixelsDistance(current_index_x, current_index_y, x, y) <toolpath_size / 2.0)
									{
										if (image_space.pixels[x][y].related_boundary_x_entry >= 0 && image_space.pixels[x][y].related_boundary_y_entry >= 0)
										{
											if (image_space.pixels[x][y].related_boundary_index_entry >= 0 && image_space.pixels[x][y].related_boundary_index_entry >= spiral_entry_affect_index - 1)
											{
												if (abs(image_space.pixels[x][y].related_boundary_index_entry - image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_index_entry) <= 1.0*toolpath_size / image_space.pixel_size)
												{
													if (image_space.pixels[x][y].bool_t&&image_space.pixels[x][y].bool_t_entry)
													{
														goon_check = false;
														break;
													}
												}
											}
										}
									}
								}
							}

							if (first_two_step_bool)
							{
								if (image_space.pixels[current_index_x][current_index_y].related_boundary_index_entry >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_index_entry >= spiral_entry_affect_index - 1)
								{
									use_this_pixel = false;
								}
							}

							if (false)
							if (spiral_exit.size()>0 && image_space.pixels[current_index_x][current_index_y].related_boundary_x_entry >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_y_entry >= 0)
							{
								if (image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_x_entry >= 0 && image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_y_entry >= 0)
								{
									if (image_space.pixels[current_index_x][current_index_y].related_boundary_index_entry >= 0 && image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_index_entry>=0)
									{
										if (image_space.pixels[current_index_x][current_index_y].related_boundary_index_entry < image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_index_entry)
										{
											use_this_pixel = false;
										}
									}
								}
							}

							if (use_this_pixel)
							{
								if (first_two_step_bool)
								{
									if (abs(p.distance_entry - toolpath_size / 2.0) < 0.0001)
									{
										perfect_index.push_back(PixelIndex(spiral_exit[spiral_exit.size() - 1].x + i, spiral_exit[spiral_exit.size() - 1].y + j));
									}

									all_index.push_back(PixelIndex(i, j));
									all_index_error.push_back(abs(p.distance_entry - toolpath_size / 2.0));
									all_index_distance.push_back(p.distance_entry);

									all_pixels.push_back(p);

									if (abs(p.distance_entry - toolpath_size / 2.0) < min_d)
									{
										min_d = abs(p.distance_entry - toolpath_size / 2.0);
										index_x = spiral_exit[spiral_exit.size() - 1].x + i;
										index_y = spiral_exit[spiral_exit.size() - 1].y + j;
										goon = true;
									}
								}
								else
								{
									if (abs(p.distance - toolpath_size / 2.0) < 0.0001)
									{
										perfect_index.push_back(PixelIndex(spiral_exit[spiral_exit.size() - 1].x + i, spiral_exit[spiral_exit.size() - 1].y + j));
									}

									all_index.push_back(PixelIndex(i, j));
									all_index_error.push_back(abs(p.distance - toolpath_size / 2.0));
									all_index_distance.push_back(p.distance);

									all_pixels.push_back(p);

									if (abs(p.distance - toolpath_size / 2.0) < min_d)
									{
										min_d = abs(p.distance - toolpath_size / 2.0);
										index_x = spiral_exit[spiral_exit.size() - 1].x + i;
										index_y = spiral_exit[spiral_exit.size() - 1].y + j;
										goon = true;
									}
								}
							}
						}
					}
				}
			}

			if (!goon_check)
			{
				goon = false;
			}

			if (goon)
			{
				if (perfect_index.size() > 1)
				{
					int spiral_index_min = (int)MAXDOUBLE;
					int spiral_index;
					for (int i = 0; i < perfect_index.size(); i++)
					{
						if (image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index_entry >= image_space.pixels[spiral_exit[spiral_exit.size() - 1].x][spiral_exit[spiral_exit.size() - 1].y].related_boundary_index_entry)
						{
							if (spiral_index_min > image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index_entry)
							{
								spiral_index_min = image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index_entry;
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

					if (spiral_exit.size() >= 2)
					{
						double max_distance_temp = image_space.PixelsDistance(index_x, index_y, spiral_exit[spiral_exit.size() - 2].x, spiral_exit[spiral_exit.size() - 2].y);

						for (int i = 0; i < all_index_error.size(); i++)
						{
							if (!(all_index[i].x == (index_x - spiral_exit[spiral_exit.size() - 1].x) && all_index[i].y == (index_y - spiral_exit[spiral_exit.size() - 1].y)))
							{
								if (abs(all_index_error[i] - min_d) < 0.0001)
								{
									double distance_temp = image_space.PixelsDistance(spiral_exit[spiral_exit.size() - 1].x + all_index[i].x, spiral_exit[spiral_exit.size() - 1].y + all_index[i].y, spiral_exit[spiral_exit.size() - 2].x, spiral_exit[spiral_exit.size() - 2].y);

									if (image_space.pixels[index_x_x][index_y_y].related_boundary_index_entry == image_space.pixels[spiral_exit[spiral_exit.size() - 1].x + all_index[i].x][spiral_exit[spiral_exit.size() - 1].y + all_index[i].y].related_boundary_index_entry)
									{
										if (max_distance_temp < distance_temp)
										{
											max_distance_temp = distance_temp;

											index_x_x = spiral_exit[spiral_exit.size() - 1].x + all_index[i].x;
											index_y_y = spiral_exit[spiral_exit.size() - 1].y + all_index[i].y;
										}
									}
									else
									{
										if (image_space.pixels[index_x_x][index_y_y].related_boundary_index_entry<image_space.pixels[spiral_exit[spiral_exit.size() - 1].x + all_index[i].x][spiral_exit[spiral_exit.size() - 1].y + all_index[i].y].related_boundary_index_entry)
										{
											max_distance_temp = distance_temp;
											index_x_x = spiral_exit[spiral_exit.size() - 1].x + all_index[i].x;
											index_y_y = spiral_exit[spiral_exit.size() - 1].y + all_index[i].y;
										}
									}
								}
							}
						}
					}

					index_x = index_x_x;
					index_y = index_y_y;
				}
			}
		

			if (goon&&!first_two_step_bool)
			{
				for (int x = index_x - toolpath_size / image_space.pixel_size; x <= index_x + toolpath_size / image_space.pixel_size; x++)
				{
					for (int y = index_y - toolpath_size / image_space.pixel_size; y <= index_y + toolpath_size / image_space.pixel_size; y++)
					{
						if (image_space.check(x, y) && image_space.PixelsDistance(index_x, index_y, x, y) <toolpath_size / 2.0 &&goon)
						{
							if (image_space.pixels[x][y].filled&&image_space.pixels[x][y].filled_entry_exit)
							{
								goon = false;
								break;
							}
						}
					}
				}
			}
			
			for (int i = 0; i < all_index.size(); i++)
			{
				if (!(index_x == spiral_exit[spiral_exit.size() - 1].x + all_index[i].x&&
					index_y == spiral_exit[spiral_exit.size() - 1].y + all_index[i].y))
				{
					image_space.pixels[spiral_exit[spiral_exit.size() - 1].x + all_index[i].x][spiral_exit[spiral_exit.size() - 1].y + all_index[i].y].failure_nb++;
				}
			}

			std::vector<PixelIndex>().swap(perfect_index);
			std::vector<PixelIndex>().swap(all_index);
			std::vector<double>().swap(all_index_error);
			std::vector<double>().swap(all_index_distance);
			std::vector<Pixel>().swap(all_pixels);
			if (goon)
			{
				InsertPixeltoSpiralExit(PixelIndex(index_x, index_y));
			}
			else
			{
				std::cout << "Over..." << std::endl;
				break;
			}

			if (!b)
			{
				break;
			}
		} while (goon);
	}

	void ToolpathGenerator::SearchOneStepforSpiralEntry(bool b)
	{
		iter_exit = 0;
		bool goon = true;

		do
		{
			//std::cerr << entry_number<<" / "<<iter_entry << std::endl;
			iter_entry++;

			double min_d = MAXDOUBLE;
			int index_x = -1;
			int index_y = -1;

			std::vector<PixelIndex> perfect_index;
			std::vector<PixelIndex> all_index;
			std::vector<double> all_index_error;
			std::vector<double> all_index_distance;
			std::vector<Pixel> all_pixels;
			bool goon_check = true;
			goon = false;

			Polygon_2 one_contour;
			bool one_contour_bool = false;
			if (spiral_entry.size() > 1)
			{
				int related_boundary_x1 = image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x_exit;
				int related_boundary_y1 = image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y_exit;

				int related_boundary_x2 = image_space.pixels[spiral_entry[spiral_entry.size() - 2].x][spiral_entry[spiral_entry.size() - 2].y].related_boundary_x_exit;
				int related_boundary_y2 = image_space.pixels[spiral_entry[spiral_entry.size() - 2].x][spiral_entry[spiral_entry.size() - 2].y].related_boundary_y_exit;

				if (related_boundary_x1 >= 0 && related_boundary_y1 >= 0 && related_boundary_x2 >= 0 && related_boundary_y2 >= 0)
				{
					int int11 = spiral_entry[spiral_entry.size() - 1].x;
					int int12 = spiral_entry[spiral_entry.size() - 1].y;

					int int21 = spiral_entry[spiral_entry.size() - 2].x;
					int int22 = spiral_entry[spiral_entry.size() - 2].y;

					one_contour.push_back(Point_2(image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].center.x, image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].center.y));
					one_contour.push_back(Point_2(image_space.pixels[spiral_entry[spiral_entry.size() - 2].x][spiral_entry[spiral_entry.size() - 2].y].center.x, image_space.pixels[spiral_entry[spiral_entry.size() - 2].x][spiral_entry[spiral_entry.size() - 2].y].center.y));

					if (related_boundary_x1 == related_boundary_x2&&related_boundary_y1 == related_boundary_y2)
					{
						one_contour.push_back(Point_2(image_space.pixels[related_boundary_x2][related_boundary_y2].center.x, image_space.pixels[related_boundary_x2][related_boundary_y2].center.y));
					}
					else
					{
						one_contour.push_back(Point_2(image_space.pixels[related_boundary_x2][related_boundary_y2].center.x, image_space.pixels[related_boundary_x2][related_boundary_y2].center.y));
						one_contour.push_back(Point_2(image_space.pixels[related_boundary_x1][related_boundary_y1].center.x, image_space.pixels[related_boundary_x1][related_boundary_y1].center.y));
					}
					one_contour_bool = true;
				}
			}

			for (int i = -1; i <= 1 && goon_check; i++)
			{
				for (int j = -1; j <= 1 && goon_check; j++)
				{
					if (image_space.check(spiral_entry[spiral_entry.size() - 1].x + i, spiral_entry[spiral_entry.size() - 1].y + j))
					{
						Pixel &p = image_space.pixels[spiral_entry[spiral_entry.size() - 1].x + i][spiral_entry[spiral_entry.size() - 1].y + j];

						int current_index_x = spiral_entry[spiral_entry.size() - 1].x + i;
						int current_index_y = spiral_entry[spiral_entry.size() - 1].y + j;
						
				

						if (p.inside&&!p.centeraxis)
						{
							bool use_this_pixel = true;

							if (image_space.pixels[current_index_x][current_index_y].failure_nb >= 2)
							{
								use_this_pixel = false;
							}

							if (spiral_entry.size() == 1 || first_two_step_bool)
							{
								for (int x = current_index_x - toolpath_size / image_space.pixel_size; x <= current_index_x + toolpath_size / image_space.pixel_size; x++)
								{
									for (int y = current_index_y - toolpath_size / image_space.pixel_size; y <= current_index_y + 1.0*toolpath_size / image_space.pixel_size; y++)
									{
										if (image_space.check(x, y) && image_space.PixelsDistance(current_index_x, current_index_y, x, y) <toolpath_size / 2.0&&use_this_pixel)
										{
											if (image_space.pixels[x][y].filled&&!image_space.pixels[x][y].filled_entry_exit)
											{
												use_this_pixel = false;
												break;
											}
										}
									}
								}
							}

							if (one_contour_bool&&one_contour.is_simple())
							{
								if (one_contour.bounded_side(Point_2(image_space.pixels[current_index_x][current_index_y].center.x, image_space.pixels[current_index_x][current_index_y].center.y)) == CGAL::ON_UNBOUNDED_SIDE)
								{
									int related_boundary_x_exit = image_space.pixels[current_index_x][current_index_y].related_boundary_x_exit;
									int related_boundary_y_exit = image_space.pixels[current_index_x][current_index_y].related_boundary_y_exit;

									int related_boundary_x1 = image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x_exit;
									int related_boundary_y1 = image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y_exit;

									if (related_boundary_x_exit == related_boundary_x1&&related_boundary_y_exit == related_boundary_y1)
									{
										Point_2 p0(image_space.pixels[current_index_x][current_index_y].center.x, image_space.pixels[current_index_x][current_index_y].center.y);
										Point_2 p1(image_space.pixels[related_boundary_x_exit][related_boundary_y_exit].center.x, image_space.pixels[related_boundary_x_exit][related_boundary_y_exit].center.y);

										Segment_2 s0(p0, p1);

										bool related_bool = false;

										std::vector<Point_2> ignore_ps;

										int ver_nb = 0;
										for (Polygon_2::Vertex_iterator ver_iter = one_contour.vertices_begin(); ver_iter != one_contour.vertices_end(); ver_iter++)
										{
											ver_nb++;
										}
										for (Polygon_2::Vertex_iterator ver_iter = one_contour.vertices_begin(); ver_iter != one_contour.vertices_end(); ver_iter++)
										{
											double  xxx1 = ver_iter->x();
											double  xxx2 = ver_iter->y();

											Polygon_2::Vertex_iterator pre_ver_iter0 = one_contour.vertices_begin();
											Polygon_2::Vertex_iterator pre_ver_iter1 = one_contour.vertices_begin() + 1;
											Polygon_2::Vertex_iterator pre_ver_iter2 = one_contour.vertices_begin() + 2;

											if (abs(ver_iter->x() - p1.x()) < 0.0000001&&abs(ver_iter->y() - p1.y()) < 0.0000001)
											{
												related_bool = true;
												ignore_ps.push_back(Point_2(p1.x(), p1.y()));

												if (ver_iter == one_contour.vertices_begin())
												{
													Polygon_2::Vertex_iterator pre_ver_iter = one_contour.vertices_begin() + ver_nb - 1;
													ignore_ps.push_back(Point_2(pre_ver_iter->x(), pre_ver_iter->y()));
												}
												else
												{
													Polygon_2::Vertex_iterator pre_ver_iter = ver_iter - 1;
													ignore_ps.push_back(Point_2(pre_ver_iter->x(), pre_ver_iter->y()));
												}

												Polygon_2::Vertex_iterator next_ver_iter = ver_iter + 1;
												if (next_ver_iter == one_contour.vertices_end())
												{
													ignore_ps.push_back(Point_2(one_contour.vertices_begin()->x(), one_contour.vertices_begin()->y()));
												}
												else
												{
													ignore_ps.push_back(Point_2(next_ver_iter->x(), next_ver_iter->y()));
												}
												break;
											}
										}

										for (Polygon_2::Edge_const_iterator edge_iter = one_contour.edges_begin(); edge_iter != one_contour.edges_end(); edge_iter++)
										{
											Segment_2 s1(edge_iter->source(), edge_iter->target());
											CGAL::Object result = CGAL::intersection(s0, s1);
											if (const Point_2 *ipoint = CGAL::object_cast<Point_2 >(&result)) {

												if (related_bool)
												{
													if (ignore_ps.size() == 3)
													{
														if (!(abs(ignore_ps[0].x() - ipoint->x()) < 0.0000001&&abs(ignore_ps[0].y() - ipoint->y()) < 0.0000001) &&
															!(abs(ignore_ps[1].x() - ipoint->x()) < 0.0000001&&abs(ignore_ps[1].y() - ipoint->y()) < 0.0000001) &&
															!(abs(ignore_ps[2].x() - ipoint->x()) < 0.0000001&&abs(ignore_ps[2].y() - ipoint->y()) < 0.0000001))
														{
															use_this_pixel = false;
															break;
														}
													}
												}
												else
												{
													use_this_pixel = false;
													break;
												}

											}
											else
												if (const Segment_2 *iseg = CGAL::object_cast<Segment_2>(&result)) {

												}
												else {

												}
										}
										std::vector<Point_2>().swap(ignore_ps);
									}
								}
								else
								{
									use_this_pixel = false;
								}
							}

							//if (false)
							if (spiral_entry.size() > 1)
							{
								int related_boundary_x_exit = image_space.pixels[current_index_x][current_index_y].related_boundary_x_exit;
								int related_boundary_y_exit = image_space.pixels[current_index_x][current_index_y].related_boundary_y_exit;

								if (image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x_exit == related_boundary_x_exit&&
									image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y_exit == related_boundary_y_exit&&
									image_space.pixels[spiral_entry[spiral_entry.size() - 2].x][spiral_entry[spiral_entry.size() - 2].y].related_boundary_x_exit == related_boundary_x_exit&&
									image_space.pixels[spiral_entry[spiral_entry.size() - 2].x][spiral_entry[spiral_entry.size() - 2].y].related_boundary_y_exit == related_boundary_y_exit)
								{
									Vector2d v0 = image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].center - image_space.pixels[related_boundary_x_exit][related_boundary_y_exit].center;
									Vector2d v1 = image_space.pixels[spiral_entry[spiral_entry.size() - 2].x][spiral_entry[spiral_entry.size() - 2].y].center - image_space.pixels[related_boundary_x_exit][related_boundary_y_exit].center;
									Vector2d v2 = image_space.pixels[current_index_x][current_index_y].center - image_space.pixels[related_boundary_x_exit][related_boundary_y_exit].center;

									if ((v1[0] * v0[1] - v1[1] * v0[0])*(v2[0] * v0[1] - v2[1] * v0[0]) > 0)
									{
										use_this_pixel = false;
									}
								}
							}

							if (spiral_entry_affect_index > 0 && !first_two_step_bool)
							{
								if (spiral_entry.size()>0 && image_space.pixels[current_index_x][current_index_y].related_boundary_x_exit >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_y_exit >= 0)
								{
									if (image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x_exit >= 0 && image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y_exit >= 0)
									{
										int int1 = image_space.pixels[image_space.pixels[current_index_x][current_index_y].related_boundary_x_exit][image_space.pixels[current_index_x][current_index_y].related_boundary_y_exit].boundary_index;
										int int2 = image_space.pixels[image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x_exit][image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y_exit].boundary_index;

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

							if (first_two_step_bool)
							for (int x = current_index_x - toolpath_size / image_space.pixel_size; x <= current_index_x + toolpath_size / image_space.pixel_size; x++)
							{
								for (int y = current_index_y - toolpath_size / image_space.pixel_size; y <= current_index_y + 1.0*toolpath_size / image_space.pixel_size; y++)
								{
									if (image_space.check(x, y) && image_space.PixelsDistance(current_index_x, current_index_y, x, y) <toolpath_size / 2.0)
									{
										if (image_space.pixels[x][y].related_boundary_x_exit >= 0 && image_space.pixels[x][y].related_boundary_y_exit >= 0)
										{
											if (image_space.pixels[x][y].related_boundary_index_exit >= 0 && image_space.pixels[x][y].related_boundary_index_exit >= spiral_exit_affect_index - 1)
											{
												Pixel &p = image_space.pixels[x][y];

												if (abs(image_space.pixels[x][y].related_boundary_index_exit - image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_index_exit)<=1.0*toolpath_size/image_space.pixel_size)
												{
													if (image_space.pixels[x][y].bool_t&&image_space.pixels[x][y].bool_t_exit)
													{
														goon_check = false;
														break;
													}
												}
											}
										}
									}
								}
							}
							
							if (first_two_step_bool)
							{
								if (image_space.pixels[current_index_x][current_index_y].related_boundary_index_exit >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_index_exit >= spiral_exit_affect_index - 1)
								{
									use_this_pixel = false;
								}
							}

							if (false)
							if (spiral_entry.size()>0 && image_space.pixels[current_index_x][current_index_y].related_boundary_x_exit >= 0 && image_space.pixels[current_index_x][current_index_y].related_boundary_y_exit >= 0)
							{
								if (image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_x_exit >= 0 && image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_y_exit >= 0)
								{
									int dsadsad = image_space.pixels[current_index_x][current_index_y].related_boundary_index_exit;
									int dsadasdas = image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_index_exit;

									if (image_space.pixels[current_index_x][current_index_y].related_boundary_index_exit >= 0 && image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_index_exit>=0)
									{
										if (image_space.pixels[current_index_x][current_index_y].related_boundary_index_exit < image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_index_exit)
										{
											use_this_pixel = false;
										}
									}
							
								}
							}

							if (use_this_pixel)
							{
								if (first_two_step_bool)
								{
									if (abs(p.distance_exit - toolpath_size / 2.0) < 0.0001)
									{
										perfect_index.push_back(PixelIndex(spiral_entry[spiral_entry.size() - 1].x + i, spiral_entry[spiral_entry.size() - 1].y + j));
									}

									all_index.push_back(PixelIndex(i, j));
									all_index_error.push_back(abs(p.distance_exit - toolpath_size / 2.0));
									all_index_distance.push_back(p.distance_exit);
									all_pixels.push_back(p);

									if (abs(p.distance_exit - toolpath_size / 2.0) < min_d)
									{
										min_d = abs(p.distance_exit - toolpath_size / 2.0);
										index_x = spiral_entry[spiral_entry.size() - 1].x + i;
										index_y = spiral_entry[spiral_entry.size() - 1].y + j;
										goon = true;
									}
								}
								else
								{
									if (abs(p.distance - toolpath_size / 2.0) < 0.0001)
									{
										perfect_index.push_back(PixelIndex(spiral_entry[spiral_entry.size() - 1].x + i, spiral_entry[spiral_entry.size() - 1].y + j));
									}

									all_index.push_back(PixelIndex(i, j));
									all_index_error.push_back(abs(p.distance - toolpath_size / 2.0));
									all_index_distance.push_back(p.distance);
									all_pixels.push_back(p);

									if (abs(p.distance - toolpath_size / 2.0) < min_d)
									{
										min_d = abs(p.distance - toolpath_size / 2.0);
										index_x = spiral_entry[spiral_entry.size() - 1].x + i;
										index_y = spiral_entry[spiral_entry.size() - 1].y + j;
										goon = true;
									}
								}
							}
						}
					}
				}
			}


			if (!goon_check)
			{
				goon = false;
			}
			
			if (goon)
			{
				if (perfect_index.size() > 1)
				{
					int spiral_index_min = (int)MAXDOUBLE;
					int spiral_index;
					for (int i = 0; i < perfect_index.size(); i++)
					{
						if (image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index_exit >= image_space.pixels[spiral_entry[spiral_entry.size() - 1].x][spiral_entry[spiral_entry.size() - 1].y].related_boundary_index_exit)
						{
							if (spiral_index_min > image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index_exit)
							{
								spiral_index_min = image_space.pixels[perfect_index[i].x][perfect_index[i].y].related_boundary_index_exit;
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

					if (spiral_entry.size() >= 2)
					{

						double max_distance_temp = image_space.PixelsDistance(index_x, index_y, spiral_entry[spiral_entry.size() - 2].x, spiral_entry[spiral_entry.size() - 2].y);

						for (int i = 0; i < all_index_error.size(); i++)
						{
							if (!(all_index[i].x == (index_x - spiral_entry[spiral_entry.size() - 1].x) && all_index[i].y == (index_y - spiral_entry[spiral_entry.size() - 1].y)))
							{
								if (abs(all_index_error[i] - min_d) < 0.0001)
								{
									double distance_temp = image_space.PixelsDistance(spiral_entry[spiral_entry.size() - 1].x + all_index[i].x, spiral_entry[spiral_entry.size() - 1].y + all_index[i].y, spiral_entry[spiral_entry.size() - 2].x, spiral_entry[spiral_entry.size() - 2].y);

									if (image_space.pixels[index_x_x][index_y_y].related_boundary_index_exit == image_space.pixels[spiral_entry[spiral_entry.size() - 1].x + all_index[i].x][spiral_entry[spiral_entry.size() - 1].y + all_index[i].y].related_boundary_index_exit)
									{
										if (max_distance_temp < distance_temp)
										{
											max_distance_temp = distance_temp;

											index_x_x = spiral_entry[spiral_entry.size() - 1].x + all_index[i].x;
											index_y_y = spiral_entry[spiral_entry.size() - 1].y + all_index[i].y;
										}
									}
									else
									{
										if (image_space.pixels[index_x_x][index_y_y].related_boundary_index_exit<image_space.pixels[spiral_entry[spiral_entry.size() - 1].x + all_index[i].x][spiral_entry[spiral_entry.size() - 1].y + all_index[i].y].related_boundary_index_exit)
										{
											max_distance_temp = distance_temp;
											index_x_x = spiral_entry[spiral_entry.size() - 1].x + all_index[i].x;
											index_y_y = spiral_entry[spiral_entry.size() - 1].y + all_index[i].y;
										}
									}
								}
							}
						}
					}

					index_x = index_x_x;
					index_y = index_y_y;
				}
			}

			if (goon&&!first_two_step_bool)
			{
				for (int x = index_x - toolpath_size / image_space.pixel_size; x <= index_x + toolpath_size / image_space.pixel_size; x++)
				{
					for (int y = index_y - toolpath_size / image_space.pixel_size; y <= index_y + toolpath_size / image_space.pixel_size; y++)
					{
						if (image_space.check(x, y) && image_space.PixelsDistance(index_x, index_y, x, y) < toolpath_size / 2.0&&goon)
						{
							if (image_space.pixels[x][y].filled&&!image_space.pixels[x][y].filled_entry_exit)
							{
								goon = false;
								break;
							}
						}
					}
				}
			}

			for (int i = 0; i < all_index.size(); i++)
			{
				if (!(index_x == spiral_entry[spiral_entry.size() - 1].x + all_index[i].x&&
					index_y == spiral_entry[spiral_entry.size() - 1].y + all_index[i].y))
				{
					image_space.pixels[spiral_entry[spiral_entry.size() - 1].x + all_index[i].x][spiral_entry[spiral_entry.size() - 1].y + all_index[i].y].failure_nb++;
				}
			}

			std::vector<PixelIndex>().swap(perfect_index);
			std::vector<PixelIndex>().swap(all_index);
			std::vector<double>().swap(all_index_error);
			std::vector<double>().swap(all_index_distance);
			std::vector<Pixel>().swap(all_pixels);
			if (goon)
			{
				InsertPixeltoSpiralEntry(PixelIndex(index_x, index_y));
			}
			else
			{
				std::cout << "Over..." << std::endl;
				break;
			}

			if (!b)
			{
				break;
			}

		} while (goon);
	}

	void ToolpathGenerator::ChooseEntryExitPoints()
	{
		image_space.Reset();
		image_space.CheckBoundary();
		image_space.UpdateDistanceFromBoundary();

		GetOnePointFromOffset(0.2, entry_point_x, entry_point_y);
		GetOnePointFromOffset(0.8, exit_point_x, exit_point_y);

		SearchUnfilledRegionBoundary(image_space.pixels[entry_point_x][entry_point_y].related_boundary_x, image_space.pixels[entry_point_x][entry_point_y].related_boundary_y);

		InsertPixeltoSpiralEntry(PixelIndex(entry_point_x, entry_point_y));
		InsertPixeltoSpiralExit(PixelIndex(exit_point_x, exit_point_y));

		SearchOneStepforSpiralEntry(true);
		SearchOneStepforSpiralExit(true);

		first_two_step_bool = true;

		for (int i = 0; i < 8; i++)
		{
		SearchOneStepforSpiralEntry(true);
		SearchOneStepforSpiralExit(true);
		}
		for (int i = 0; i < input_int_2; i++)
			SearchOneStepforSpiralEntry(false);
		return;

		do
		{
			std::cerr << "Whole iter: " << entry_number << std::endl;

			entry_number++;

			int spiral_entry_nb = spiral_entry.size(); 
			int spiral_exit_nb = spiral_exit.size();
			SearchOneStepforSpiralEntry(true);
			SearchOneStepforSpiralExit(true);

			if (spiral_entry_nb == spiral_entry.size() && spiral_exit_nb == spiral_exit.size())
			break;

		} while (true);

		std::vector<PixelIndex> insertpixels;
		ConnectTwoPixels(spiral_entry[spiral_entry.size() - 1], spiral_exit[spiral_exit.size() - 1], insertpixels);
		
		for (int i = 0; i < insertpixels.size(); i++)
		{
			spiral_entry.push_back(insertpixels[i]);
		}
	}

}




