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
		draw_aixs =false;
		draw_contour = true;
		draw_offsets = false;
		draw_entry_exit_spiral = true;
		draw_turning_points = true;
		draw_spiral = true;

		line_width = 5;
		point_size = 5;

		entry_d_0 = 0.8;
		exit_d_0 = 0.2;

		contour_size = 30.0; //mm
		
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
		file >> load_path;

		file >> str;
		file >> draw_spiral;

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
		file >> smooth_number;

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
		}
		
		if (work_model == 2)
		{

			ArchinedeanSpiral(contour);
		}

		if (work_model == 3)
		{

			GenerateZigzag();
		}

		if (work_model == 4)
		{
			ArchinedeanSpiralTrick(contour);
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

			OutputPath(one_single_path, output_str + ".dat");

			Output_Obj_1(output_str + ".obj");
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
				one_contour.push_back(Point_2(x,y));
			}
			one_contour.reverse_orientation();
			if (i == 0)
				poly_contours = Polygon_with_holes(one_contour);
			else
				poly_contours.add_hole(one_contour);

		}
		file.clear();
		file.close();

		double center_x = (poly_contours.bbox().xmax() + poly_contours.bbox().xmin())/2.0;
		double center_y = (poly_contours.bbox().ymax() + poly_contours.bbox().ymin())/2.0;

		double scale = (poly_contours.bbox().ymax() - poly_contours.bbox().ymin())>(poly_contours.bbox().xmax() - poly_contours.bbox().xmin()) ? 
			(poly_contours.bbox().ymax() - poly_contours.bbox().ymin()) : (poly_contours.bbox().xmax() - poly_contours.bbox().xmin());
		scale = scale /4.5;

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
			x += 1.0;
			y += 0.1;
			one_contour.push_back(Point_2(x,y));
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
				x += 1.0;
				y += 0.1;
				one_contour.push_back(Point_2(x, y));
			}
			contours.add_hole(one_contour);
		}
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
			glLineWidth(line_width);
			glColor3f(65.0 / 255.0, 128.0 / 255.0, 187.0 / 255.0);

			glBegin(GL_LINE_LOOP);
			for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
			{
				glVertex3f(ver_iter->x(), ver_iter->y(), 0.0);
			}
			glEnd();

			for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
			{
				glBegin(GL_LINE_LOOP);
				for (Polygon_2::Vertex_iterator ver_iter = hole_iter->vertices_begin(); ver_iter != hole_iter->vertices_end(); ver_iter++)
				{
					glVertex3f(ver_iter->x(), ver_iter->y(), 0.0);
				}
				glEnd();
			}
		}

		//draw offset
		if (draw_offsets)
		{
			glLineWidth(1);
			//glColor3f(0.5, 0.5, 0.0);
			glColor3f(0.75, 0.75, 0.75);

			for (int i = 0; i < offsets.size(); i++)
			{
				glBegin(GL_LINE_LOOP);
				for (int j = 0; j < offsets[i].size(); j++)
				{
					glVertex3f(offsets[i][j][0], offsets[i][j][1], 0.0);
				}
				glEnd();
			}


			if (debug_int_0 == 100)
			{
				glPointSize(point_size);
				glColor3f(0.0, 0.0, 1.0);
				glBegin(GL_POINTS);
				for (int i = 0; i < offsets.size(); i++)
				{
					for (int j = 0; j < offsets[i].size(); j++)
					{
						glVertex3f(offsets[i][j][0], offsets[i][j][1], 0.0);
					}
				}
				glEnd();
			}
			
			if (false)
			for (int i = 0; i < decompose_offset.size(); i++)
			{
				double color_0 = Random();
				double color_1 = Random();
				double color_2 = Random();

				for (int j = 0; j < decompose_offset[i].size(); j++)
				{
					glLineWidth(6);
					//glColor3f(0.2, 0.2, 0.0);
					glColor3f(color_0, color_1, color_2);

					glBegin(GL_LINE_LOOP);
					for (int m = 0; m < offsets[decompose_offset[i][j]].size(); m++)
					{
						glVertex3f(offsets[decompose_offset[i][j]][m][0], offsets[decompose_offset[i][j]][m][1], 0.0);
					}
					glEnd();
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

		if (draw_entry_exit_spiral)
		{
			if (work_model == 2 || work_model == 3 || work_model == 4 || work_model == 5 || work_model == 6)
			{
				glLineWidth(line_width);
				glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
				glBegin(GL_LINE_STRIP);
				for (int j = 0; j < entry_spiral.size(); j++)
				{
					glVertex3f(entry_spiral[j][0], entry_spiral[j][1], 0.0);
				}
				glEnd();


				glLineWidth(line_width);
				glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
				glBegin(GL_LINE_STRIP);
				for (int j = 0; j < exit_spiral.size(); j++)
				{
					glVertex3f(exit_spiral[j][0], exit_spiral[j][1], 0.0);
				}
				glEnd();
			}


			if (work_model == 7)
			{
				if (draw_spiral)
				{
					for (int i = 0; i < entry_spirals.size(); i++)
					{
						glLineWidth(line_width);
						glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
						glBegin(GL_LINE_STRIP);
						for (int j = 0; j < entry_spirals[i].size(); j++)
						{
							glVertex3f(entry_spirals[i][j][0], entry_spirals[i][j][1], 0.0);
						}
						glEnd();
					}

					for (int i = 0; i < exit_spirals.size(); i++)
					{
						glLineWidth(line_width);
						glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
						glBegin(GL_LINE_STRIP);
						for (int j = 0; j < exit_spirals[i].size(); j++)
						{
							glVertex3f(exit_spirals[i][j][0], exit_spirals[i][j][1], 0.0);
						}
						glEnd();
					}
				}
				
				//pathes
				if (false)
				{
					for (int i = 0; i < pathes.size(); i++)
					{

						glLineWidth(line_width);
						glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
						glBegin(GL_LINE_STRIP);
						for (int j = 0; j < pathes[i].size(); j++)
						{
							glVertex3f(pathes[i][j][0], pathes[i][j][1], 0.0);
						}
						glEnd();
					}
				}
			}
		}


		//one_single_path
		if (iiiiii % 2 == 1)
		{
			glLineWidth(line_width*1.2);
			glColor3f(0.0, 0.0, 1.0);
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j <one_single_path.size(); j++)
			{
				glVertex3f(one_single_path[j][0], one_single_path[j][1], 0.0);
			}
			glEnd();
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


		glPointSize(point_size*2);
		glColor3f(0.0, 1.0, 1.0);
		glBegin(GL_POINTS);
		for (int j = 0; j < debug_points.size(); j++)
		{
			glVertex3f(debug_points[j][0], debug_points[j][1], 0.0);
		}
		glEnd();

		glLineWidth(line_width*2);
		glColor3f(0.0, 0.0, 1.0);

		
		for (int i = 0; i < debug_lines.size(); i++)
		{
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j <debug_lines[i].size(); j++)
			{
				glVertex3f(debug_lines[i][j][0], debug_lines[i][j][1], 0.0);
			}
			glEnd();
		}

	
		if (draw_turning_points)
		{
			glPointSize(point_size);
			glColor3f(0.0,1.0,0.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < turning_points_exit.size(); i++)
			{
				glVertex3f(turning_points_exit[i][0], turning_points_exit[i][1], 0.0);
			}
			glEnd();

			glPointSize(point_size);
			glColor3f(1.0,0.0,0.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < turning_points_entry.size(); i++)
			{
				glVertex3f(turning_points_entry[i][0], turning_points_entry[i][1], 0.0);
			}
			glEnd();


		}



	}


}