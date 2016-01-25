#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <iostream>
#include <fstream>
#include <list>
#include <set>

#include "Circuit.h"

//#include "ToolPathTimeEstimator.hpp"

namespace hpcg {

	ToolpathGenerator::ToolpathGenerator(){

	}

	int iiiiii = 0;
	void  ToolpathGenerator::StepDebug()
	{
		iiiiii++;
	}

	double Random()
	{
		return ((double)(rand() % 1000)) / 1000.0;
	}


	void ToolpathGenerator::Init(TMeshProcessor* render)
	{
		srand(0);

		image_space.pixel_size = 0.05;
		image_space.toolpath_size = 0.2;
		toolpath_size = 0.2;

		draw_pixels = false;
		draw_aixs = false;
		draw_contour = true;
		draw_offsets = false;
		draw_entry_exit_spiral = true;
		draw_turning_points = true;
		draw_spiral = true;

		sampling_one_single_path = false;

		graph_scale = 4.5;
		graph_move_x = 1.0;
		graph_move_y = 0.1;

		line_width = 5;
		point_size = 5;

		entry_d_0 = 0.8;
		exit_d_0 = 0.2;

		contour_size = 30.0; //mm


		contour_path_size = 3;
		contour_color_r = 1.0;
		contour_color_g = 1.0;
		contour_color_b = 1.0;

		std::string input_path;

		std::ifstream inputfile("D:\\test.txt", std::ios::in);

		int input_int;
		inputfile >> input_int;

		for (int i = 0; i < input_int; i++)
		{
			inputfile >> input_path;
		}

		inputfile.clear();
		inputfile.close();

		std::string path = input_path;
		std::ifstream file(path, std::ios::in);

		if (!file)
		{
			std::cout << "" << std::endl;
			return;
		}
		std::string str;
		std::string output_str;

		file >> debug_int_0;
		file >> debug_int_1;
		file >> debug_int_2;

		file >> str;
		file >> graph_scale;

		file >> str;
		file >> graph_move_x;

		file >> str;
		file >> graph_move_y;

		file >> str;
		file >> load_path;

		file >> str;
		file >> draw_spiral;

		file >> str;
		file >> draw_contour;

		file >> str;
		file >> draw_turning_points;

		file >> str;
		file >> draw_offsets;

		file >> str;
		file >> line_width;

		file >> str;
		file >> point_size;

		file >> str;
		file >> work_model;

		file >> str;
		file >> sampling_one_single_path;

		file >> str;
		file >> smooth_number;

		//smooth_boundary
		file >> str;
		file >> smooth_boundary;

		file >> str;
		file >> contour_size;

		file >> str;
		file >> input_toolpath_size;

		file >> str;
		file >> entry_d_0;

		file >> str;
		file >> exit_d_0;

		file >> str;
		file >> output_str;

		file >> str;
		file >> use_save_offset_file;

		file >> str;
		file >> offset_file;


		file >> str;
		file >> peter_data_path;

		
		//int contour_path_size;
		//double contour_color_r;
		//double contour_color_g;
		//double contour_color_b;

		//int filling_path_size;
		//double filling_color_r;
		//double filling_color_g;
		//double filling_color_b;

		file >> str;
		file >> contour_path_size;
		file >> str;
		file >> contour_color_r;
		file >> str;
		file >> contour_color_g;
		file >> str;
		file >> contour_color_b;

		file >> str;
		file >> filling_path_size;
		file >> str;
		file >> filling_color_r;
		file >> str;
		file >> filling_color_g;
		file >> str;
		file >> filling_color_b;


		/*
		double entry_spiral_size;
		double entry_spiral_color_r;
		double entry_spiral_color_g;
		double entry_spiral_color_b;

		double exit_spiral_size;
		double exit_spiral_color_r;
		double exit_spiral_color_g;
		double exit_spiral_color_b;
		*/
		file >> str;
		file >> entry_spiral_size;
		file >> str;
		file >> entry_spiral_color_r;
		file >> str;
		file >> entry_spiral_color_g;
		file >> str;
		file >> entry_spiral_color_b;

		file >> str;
		file >> exit_spiral_size;
		file >> str;
		file >> exit_spiral_color_r;
		file >> str;
		file >> exit_spiral_color_g;
		file >> str;
		file >> exit_spiral_color_b;


		file.clear();
		file.close();

		OutputPathTwoCircles();

		m_render = render;
		LoadContour();

		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		double max_double;

		if (abs(contours.bbox().xmax() - contours.bbox().xmin()) >
			abs(contours.bbox().ymax() - contours.bbox().ymin()))
		{
			max_double = abs(contours.bbox().xmax() - contours.bbox().xmin());
		}
		else
		{
			max_double = abs(contours.bbox().ymax() - contours.bbox().ymin());
		}

		toolpath_size = input_toolpath_size / contour_size*max_double;

		image_space.toolpath_size = toolpath_size;

		//FillingAlgorithm();
		//FermatSpiral();

		if (work_model == 1)
		{
			for (int i = 0; i < smooth_number; i++)
			{
				PolygonSmoothing();
				DirectlyPolygonSmoothing();
			}
		}

		if (work_model == 2)
		{

			FermatSpiral(contour, entry_d_0, exit_d_0);
			//ArchinedeanSpiral(contour);

			if (sampling_one_single_path)
			{
				std::vector<Vector2d>().swap(one_single_path);

				for (int i = 0; i < entry_spiral.size(); i++)
				{
					one_single_path.push_back(entry_spiral[i]);
				}

				for (int i = exit_spiral.size() - 1; i >= 0; i--)
				{
					one_single_path.push_back(exit_spiral[i]);
				}

				std::vector<Vector2d> temp;
				for (int i = 0; i < one_single_path.size(); i++)
				{
					temp.push_back(one_single_path[i]);
				}

				std::vector<Vector2d>().swap(one_single_path);
				for (int i = 0; i < 2000; i++)
				{
					one_single_path.push_back(Strip::GetOnePointFromStrip((double)i / 2000.0, temp));
				}

				OutputPath(one_single_path, output_str + ".dat");
				Output_Obj_1(output_str + ".obj");
			}
		}

		if (work_model == 3)
		{
			GenerateZigzag();
		}

		if (work_model == 4)
		{
			ArchinedeanSpiral(contour);
		}

		if (work_model == 5)
		{
			OutputPath_Hollow(output_str);
		}

		if (work_model == 6)
		{
			TurnPath();
			return;
		}

		if (work_model == 7)
		{
			FillingAlgorithmBasedOnOffsets();

			if (sampling_one_single_path)
				if (one_single_path.size() != 0)
				{
					std::vector<Vector2d> temp;
					for (int i = 0; i < one_single_path.size(); i++)
					{
						temp.push_back(one_single_path[i]);
					}

					std::vector<Vector2d>().swap(one_single_path);
					for (int i = 0; i <= 10000; i++)
					{
						one_single_path.push_back(Strip::GetOnePointFromStrip((double)i / 10000.0, temp));
					}

					BoundaryTag();

					OutputPath(one_single_path, output_str + ".dat");
					Output_Obj_1(output_str + ".obj");
				}
		}


		if (work_model == 8)
		{
			for (int i = 0; i < smooth_number; i++)
			{
				PolygonSmoothing();
				DirectlyPolygonSmoothing();
			}

			InputOneSinglePath(peter_data_path);
		}

		//ComputeOffsets_temp();
		//OutputPath(one_single_path, output_str);
		//OutputPath(entry_spiral, str);

		std::vector<Vector2d>().swap(contour);
	}

	void ToolpathGenerator::LoadContour()
	{
		//load
		std::ifstream file(load_path, std::ios::in);

		if (!file)
		{
			std::cout << "" << std::endl;
			return;
		}

		int contour_number;
		file >> contour_number;

		Polygon_with_holes poly_contours;
		for (int i = 0; i < contour_number; i++){
			int side_number;
			file >> side_number;

			Polygon_2 one_contour;
			for (int j = 0; j < side_number; j++){
				double x, y;
				file >> x >> y;
				one_contour.push_back(Point_2(x, y));
			}
			one_contour.reverse_orientation();
			if (i == 0)
				poly_contours = Polygon_with_holes(one_contour);
			else
				poly_contours.add_hole(one_contour);

		}
		file.clear();
		file.close();

		double center_x = (poly_contours.bbox().xmax() + poly_contours.bbox().xmin()) / 2.0;
		double center_y = (poly_contours.bbox().ymax() + poly_contours.bbox().ymin()) / 2.0;

		double scale = (poly_contours.bbox().ymax() - poly_contours.bbox().ymin())>(poly_contours.bbox().xmax() - poly_contours.bbox().xmin()) ?
			(poly_contours.bbox().ymax() - poly_contours.bbox().ymin()) : (poly_contours.bbox().xmax() - poly_contours.bbox().xmin());
		scale = scale / graph_scale;

		Polygon_2 one_contour;
		for (Polygon_2::Vertex_iterator ver_iter = poly_contours.outer_boundary().vertices_begin(); ver_iter != poly_contours.outer_boundary().vertices_end(); ver_iter++)
		{
			double x = ver_iter->x();
			double y = ver_iter->y();
			x -= center_x;
			y -= center_y;
			x = -x;
			y = -y;
			x = x / scale;
			y = y / scale;
			x += graph_move_x;
			y += graph_move_y;
			one_contour.push_back(Point_2(x, y));
		}

		if (!one_contour.is_clockwise_oriented())
		{
			one_contour.reverse_orientation();
		}


		contours = Polygon_with_holes(one_contour);

		for (Polygon_with_holes::Hole_iterator hole_iter = poly_contours.holes_begin(); hole_iter != poly_contours.holes_end(); hole_iter++)
		{
			Polygon_2 one_contour;
			for (Polygon_2::Vertex_iterator ver_iter = hole_iter->vertices_begin(); ver_iter != hole_iter->vertices_end(); ver_iter++)
			{
				double x = ver_iter->x();
				double y = ver_iter->y();
				x -= center_x;
				y -= center_y;
				x = -x;
				y = -y;
				x = x / scale;
				y = y / scale;
				x += graph_move_x;
				y += graph_move_y;
				one_contour.push_back(Point_2(x, y));
			}

			if (one_contour.is_clockwise_oriented())
			{
				one_contour.reverse_orientation();
			}
			contours.add_hole(one_contour);
		}
	}


	void ToolpathGenerator::Draw_A_Line(Vector2d v0, Vector2d v1, double width, double color_r, double color_g, double color_b)
	{
		if (true)
		//if (debug_int_0 == 0)
		{
			glLineWidth(line_width);
			glColor3f(color_r, color_g, color_b);
			glBegin(GL_LINE);
			glVertex3f(v0[0], v0[1], 0.0);
			glVertex3f(v1[0], v1[1], 0.0);
			glEnd();


			glLineWidth(line_width * 2);
			glColor3f(color_r, color_g, color_b);
			glBegin(GL_LINES);
			glVertex3f(v0[0], v0[1], 0.0);
			glVertex3f(v1[0], v1[1], 0.0);

			glEnd();

			return;
		}

		Vector2d n(v0[0] - v1[0], v0[1] - v1[1]);
	
		if (debug_int_0<0)
		{
			double length0 = sqrt(n[0] * n[0] + n[1] * n[1]);
			n[0] = n[0] / length0;
			n[1] = n[1] / length0;

			v0[0] += n[0] * width / 2.0;
			v0[1] += n[1] * width / 2.0;

			v1[0] -= n[0] * width / 2.0;
			v1[1] -= n[1] * width / 2.0;
		}


		Vector2d pn(-n[1],n[0]);
		double length = sqrt(pn[0] * pn[0] + pn[1] * pn[1]);
		pn[0] = pn[0] / length;
		pn[1] = pn[1] / length;

		Vector2d p0(v0[0] + width / 2.0*pn[0], v0[1] + width / 2.0*pn[1]);
		Vector2d p1(v1[0] + width / 2.0*pn[0], v1[1] + width / 2.0*pn[1]);
		Vector2d p2(v0[0] - width / 2.0*pn[0], v0[1] - width / 2.0*pn[1]);
		Vector2d p3(v1[0] - width / 2.0*pn[0], v1[1] - width / 2.0*pn[1]);

		glColor3f(color_r, color_g, color_b);

		glBegin(GL_POLYGON);

		glVertex3f(p0[0], p0[1], 0.0);
		glVertex3f(p1[0], p1[1], 0.0);
		glVertex3f(p3[0], p3[1], 0.0);
		glVertex3f(p2[0], p2[1], 0.0);

		glEnd();

	}

	void ToolpathGenerator::Rendering()
	{


		//draw pixels
		if (draw_pixels)
		{
			glLineWidth(1);
			glColor3f(0.1, 0.1, 0.1);

			for (int i = 0; i < image_space.pixels.size(); i++)
			{
				for (int j = 0; j < image_space.pixels[i].size(); j++)
				{
					if (image_space.pixels[i][j].inside)
					{
						if (image_space.pixels[i][j].filled)
						{
							glColor3f(0.8, 0.8, 0.8);
							glBegin(GL_POLYGON);
							DrawPixel(i, j);
							glEnd();
						}
					}
				}
			}
		}


		//draw contour
		if (draw_contour)
		{

			std::vector<Vector2d> one_path;

			for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
			{
				one_path.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
			}

			for (int i = 0; i < one_path.size(); i++)
			{
				Draw_A_Line(one_path[i], one_path[(i + 1) % one_path.size()], contour_path_size, contour_color_r, contour_color_g, contour_color_b);
			}

			for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
			{

				std::vector<Vector2d>().swap(one_path);

				for (Polygon_2::Vertex_iterator ver_iter = hole_iter->vertices_begin(); ver_iter != hole_iter->vertices_end(); ver_iter++)
				{
					one_path.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				}
				for (int i = 0; i < one_path.size(); i++)
				{
					Draw_A_Line(one_path[i], one_path[(i + 1) % one_path.size()], contour_path_size, contour_color_r, contour_color_g, contour_color_b);
				}
			}

			//Draw_A_Line(Vector2d(0.0, 0.0), Vector2d(0.5, 0.0), 0.1, 1.0, 0.0, 0.0);
		}


		//draw offset
		if (draw_offsets)
		{
			for (int i = 0; i < offsets.size(); i++)
			{
				for (int j = 0; j < offsets[i].size(); j++)
				{
					Draw_A_Line(offsets[i][j], offsets[i][(j + 1) % offsets[i].size()], filling_path_size, 0.9, 0.9, 0.9);
				}
			}

			if (draw_turning_points)
			for (int i = 0; i < decompose_offset.size(); i++)
			{
				double color_0 = Random();
				double color_1 = Random();
				double color_2 = Random();

				for (int j = 0; j < decompose_offset[i].size(); j++)
				{
					//glColor3f(color_0, color_1, color_2);

					for (int m = 0; m < offsets[decompose_offset[i][j]].size(); m++)
					{
						Draw_A_Line(offsets[decompose_offset[i][j]][m], offsets[decompose_offset[i][j]][(m + 1) % offsets[decompose_offset[i][j]].size()], filling_path_size, color_0, color_1, color_2);
					}

				}
			}
		}

		//draw axis
		if (draw_aixs)
		{
			glLineWidth(1);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_LINES);
			glVertex3f(-1000.0, 0.0, 0.0);
			glVertex3f(1000.0, 0.0, 0.0);
			glEnd();

			glColor3f(0.0, 0.0, 1.0);
			glBegin(GL_LINES);
			glVertex3f(0.0, -1000.0, 0.0);
			glVertex3f(0.0, 1000.0, 0.0);
			glEnd();
		}

		if (draw_spiral)
		{
			if (iiiiii % 2 == 0)
			if (work_model == 2 || work_model == 3 || work_model == 4 || work_model == 5 || work_model == 6)
			{
				if (entry_spiral.size() != 0)
				for (int j = 0; j < entry_spiral.size()-1; j++)
				{
					Draw_A_Line(entry_spiral[j], entry_spiral[(j + 1) % entry_spiral.size()], entry_spiral_size, entry_spiral_color_r, entry_spiral_color_g, entry_spiral_color_b);
				}
			
				if (exit_spiral.size()!=0)
				for (int j = 0; j < exit_spiral.size()-1; j++)
				{
					Draw_A_Line(exit_spiral[j], exit_spiral[(j + 1) % exit_spiral.size()], exit_spiral_size, exit_spiral_color_r, exit_spiral_color_g, exit_spiral_color_b);
				}
			}

			if (work_model == 7)
			{
				if (iiiiii % 2 == 0)
				{
					for (int i = 0; i < entry_spirals.size(); i++)
					{
						if (entry_spirals[i].size()!=0)
						for (int j = 0; j < entry_spirals[i].size()-1; j++)
						{
							Draw_A_Line(entry_spirals[i][j], entry_spirals[i][(j + 1) % entry_spirals[i].size()], entry_spiral_size, entry_spiral_color_r, entry_spiral_color_g, entry_spiral_color_b);
						}
					}
					for (int i = 0; i < exit_spirals.size(); i++)
					{
						if (exit_spirals[i].size() != 0)
						for (int j = 0; j < exit_spirals[i].size()-1; j++)
						{
							Draw_A_Line(exit_spirals[i][j], exit_spirals[i][(j + 1) % exit_spirals[i].size()], exit_spiral_size, exit_spiral_color_r, exit_spiral_color_g, exit_spiral_color_b);
						}
					}
				}
			}

			//one_single_path
			if (iiiiii % 2 == 1)
			{
				if (debug_int_0 == -2 && debug_int_2 == -2)
				{
					glPointSize(point_size);
					glColor3f(0.0, 0.0, 1.0);
					glBegin(GL_POINTS);
					for (int j = 0; j <one_single_path.size(); j++)
					{
						glVertex3f(one_single_path[j][0], one_single_path[j][1], 0.0);


					}
					glEnd();
				}
				else
				{
					if (one_single_path.size() != 0)
						for (int j = 0; j <one_single_path.size() - 1; j++)
						{
							Draw_A_Line(one_single_path[j], one_single_path[(j + 1) % one_single_path.size()], filling_path_size*2.0, filling_color_r, filling_color_g, filling_color_b);

						}

					if (false)
						if (one_single_path_boundary.size() != 0)
						{
							glPointSize(point_size);
							glColor3f(1.0, 0.0, 0.0);
							glBegin(GL_POINTS);
							for (int j = 0; j <one_single_path.size(); j++)
							{
								if (one_single_path_boundary[j])
									glVertex3f(one_single_path[j][0], one_single_path[j][1], 0.0);
							}
							glEnd();
						}
				}
			}
		}



		//pathes_temp
		if (iiiiii % 2 == 0)
		{
			for (int i = 0; i < pathes_temp.size(); i++)
			{
				if (pathes_temp[i].size() == 1)
				{
					glPointSize(point_size);
					glColor3f(0.0, 0.0, 1.0);
					glBegin(GL_POINTS);
					for (int j = 0; j < pathes_temp[i].size(); j++)
					{
						glVertex3f(pathes_temp[i][j][0], pathes_temp[i][j][1], 0.0);
					}
					glEnd();
				}
				else
				{
					glLineWidth(line_width);
					glColor3f(0.0, 0.0, 1.0);
					glBegin(GL_LINE_STRIP);
					for (int j = 0; j < pathes_temp[i].size(); j++)
					{
						glVertex3f(pathes_temp[i][j][0], pathes_temp[i][j][1], 0.0);
					}
					glEnd();
				}
			}
		}

		glLineWidth(line_width * 2);
		glColor3f(0.0, 1.0, 1.0);
		for (int i = 0; i < debug_lines.size(); i++)
		{
			glBegin(GL_LINES);
			for (int j = 0; j <debug_lines[i].size(); j++)
			{
				glVertex3f(debug_lines[i][j][0], debug_lines[i][j][1], 0.0);
			}
			glEnd();
		}

		glPointSize(point_size * 2);
		glColor3f(0.0, 1.0, 1.0);
		glBegin(GL_POINTS);
		for (int j = 0; j < debug_points.size(); j++)
		{
			glVertex3f(debug_points[j][0], debug_points[j][1], 0.0);
		}
		glEnd();


	}


}