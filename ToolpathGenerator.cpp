#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <iostream>
#include <fstream>
#include <list>
#include <set>


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
		draw_voronoi = false;
		draw_medial_axis = false;
		draw_minimal_points = false;
		draw_maximal_points = false;
		draw_entry_exit_points = true;
		draw_cutting_points = false;

		draw_polygons_entry_exit = false;
		draw_inner_concave_points = false;

		draw_save_critical_points = false;

		line_width = 5;
		point_size = 5;

		entry_d_0 = 0.8;
		exit_d_0 = 0.2;

		std::string path = "D:\\test.txt";
		std::ifstream file(path, std::ios::in);

		if (!file)
		{
			std::cout << "" << std::endl;
			return;
		}
		std::string str;

		file >> debug_int_0;
		file >> debug_int_1;

		file >> str;
		file >> draw_medial_axis;

		file >> str;
		file >> draw_minimal_points;

		file >> str;
		file >> draw_maximal_points;

		file >> str;
		file >> load_path;

		file >> str;
		file >> draw_entry_exit_spiral;

		file >> str;
		file >> draw_turning_points;


		file >> str;
		file >> draw_offsets;

		file >> str;
		file >> draw_cutting_points;


		file >> str;
		file >> draw_entry_exit_points;

		file >> str;
		file >> line_width;

		file >> str;
		file >> point_size;

		file >> str;
		file >> work_model;

		file >> str;
		file >> smooth_number;

		file >> str;
		file >> toolpath_size;

		file >> str;
		file >> entry_d_0;

		file >> str;
		file >> exit_d_0;

		file >> str;

		file.clear();
		file.close();
		
		image_space.toolpath_size = toolpath_size;

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

	
		//FillingAlgorithm();
		//FermatSpiral();

		if (work_model == 1)
		{
			FillingAlgorithm();
		}
		
		if (work_model == 2)
		{

			toolpath_size = 0.4 / 30.0*max_double;
			ArchinedeanSpiral(contour);
		}

		if (work_model == 3)
		{

			toolpath_size = 0.4 / 30.0*max_double;
			GenerateZigzag();
		}

		if (work_model == 4)
		{
			//toolpath_size = 0.4 / 30.0*max_double;
			toolpath_size =0.4 / 30.0*max_double;
			
			ArchinedeanSpiralTrick(contour);
		}

		if (work_model == 5)
		{
			toolpath_size = 0.4 / 30.0*max_double;
			GenerateZigzagForCircle();
		}

		if (work_model == 6)
		{
			toolpath_size = 0.4 / 30.0*max_double;
			ArchinedeanSpiralTrickForCircle(contour);
		}

		if (work_model == 7)
		{
			FillingAlgorithmBasedOnOffsets();
		}

		//ComputeOffsets_temp();
		  
		OutputPath(entry_spiral, str);

		//OutputPath(str);
		
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
				if (false)
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
				if (true)
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

				//pathes_temp
				if (iiiiii%2==1)
					for (int i = 0; i < pathes_temp.size(); i++)
					{

						glLineWidth(line_width);
						glColor3f(1.0,0.0,0.0);
						glBegin(GL_LINE_STRIP);
						for (int j = 0; j < pathes_temp[i].size(); j++)
						{
							glVertex3f(pathes_temp[i][j][0], pathes_temp[i][j][1], 0.0);
						}
						glEnd();
					}

			}

	
			for (int i = 0; i < region.entry_spirals.size(); i++)
			{
				glLineWidth(line_width);
				glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
				glBegin(GL_LINE_STRIP);
				for (int j = 0; j < region.entry_spirals[i].size(); j++)
				{
					glVertex3f(region.entry_spirals[i][j][0], region.entry_spirals[i][j][1], 0.0);
				}
				glEnd();

				glLineWidth(line_width);
				glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
				glBegin(GL_LINE_STRIP);
				for (int j = 0; j < region.exit_spirals[i].size(); j++)
				{
					glVertex3f(region.exit_spirals[i][j][0], region.exit_spirals[i][j][1], 0.0);
				}
				glEnd();

				if (false)
				{
					glLineWidth(line_width);
					glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
					glBegin(GL_LINE_STRIP);

					if (region.entry_spirals[i].size() > 0 && region.exit_spirals[i].size())
					{
						glVertex3f(region.entry_spirals[i][region.entry_spirals[i].size() - 1][0], region.entry_spirals[i][region.entry_spirals[i].size() - 1][1], 0.0);
						glVertex3f(region.exit_spirals[i][region.exit_spirals[i].size() - 1][0], region.exit_spirals[i][region.exit_spirals[i].size() - 1][1], 0.0);
					}
					glEnd();
				}
			}

			if (false)
			for (int i = 0; i < region.connected_regions.size(); i = i + 2)
			{
				glLineWidth(line_width);
				glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
				glBegin(GL_LINE_STRIP);
				glVertex3f(region.connected_regions[i][0], region.connected_regions[i][1], 0.0);
				glVertex3f(region.connected_regions[i+1][0], region.connected_regions[i+1][1], 0.0);
				glEnd();
			}

		}


		if (draw_entry_exit_points)
		{

			for (int i = 0; i < region.polygons_entry_exit.size(); i++)
			{
				for (int j = 0; j < region.polygons_entry_exit[i].size(); j++)
				{
					glPointSize(point_size);
					glColor3f(200 / 255.0, 78 / 255.0,46.0 / 255.0);
					glBegin(GL_POINTS);
					glVertex3f(region.polygons_entry_exit[i][j][0], region.polygons_entry_exit[i][j][1], 0.0);
					glEnd();
				}
			}

			if (false)
			if (entry_spiral.size() > 0 && exit_spiral.size())
			{
				glPointSize(point_size);
				glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
				glBegin(GL_POINTS);
				glVertex3f(entry_spiral[0][0], entry_spiral[0][1], 0.0);
				glEnd();

				glPointSize(point_size);
				glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
				glBegin(GL_POINTS);
				glVertex3f(exit_spiral[0][0], exit_spiral[0][1], 0.0);
				glEnd();
			}
		}

		if (draw_voronoi)
		{
			for (int i = 0; i < region.sdg.voronoi_edge_points.size(); i = i + 2)
			{
				glLineWidth(2);
				glColor3f(0.5, 0.5, 0.0);
				glBegin(GL_LINE_STRIP);
				glVertex3f(region.sdg.voronoi_edge_points[i][0], region.sdg.voronoi_edge_points[i][1], 0.0);
				glVertex3f(region.sdg.voronoi_edge_points[(i + 1) % region.sdg.voronoi_edge_points.size()][0], region.sdg.voronoi_edge_points[(i + 1) % region.sdg.voronoi_edge_points.size()][1], 0.0);
				glEnd();
			}

		}

		if (draw_medial_axis)
		{
			
			for (int i = 0; i < region.sdg.medial_axis_points.size(); i = i + 2)
			{
				glLineWidth(2 * 2);
				glColor3f(0.4, 0.4, 0.0);
				glBegin(GL_LINE_STRIP);
				glVertex3f(region.sdg.medial_axis_points[i][0], region.sdg.medial_axis_points[i][1], 0.0);
				glVertex3f(region.sdg.medial_axis_points[(i + 1) % region.sdg.medial_axis_points.size()][0], region.sdg.medial_axis_points[(i + 1) % region.sdg.medial_axis_points.size()][1], 0.0);
				glEnd();
			}
			if (false)
				for (int i = 0; i < region.sdg.medial_axis_points.size(); i = i + 2)
				{
					glPointSize(2 * 2 + 2);
					glColor3f(1.0, 0.0, 0.0);
					glBegin(GL_POINTS);
					glVertex3f(region.sdg.medial_axis_points[i][0], region.sdg.medial_axis_points[i][1], 0.0);
					glVertex3f(region.sdg.medial_axis_points[(i + 1) % region.sdg.medial_axis_points.size()][0], region.sdg.medial_axis_points[(i + 1) % region.sdg.medial_axis_points.size()][1], 0.0);
					glEnd();
				}
		}

		if (draw_cutting_points)
		{
			for (int i = 0; i < region.cutting_points.size(); i = i + 2)
			{
				glPointSize(point_size);
				glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
				glBegin(GL_LINE_STRIP);
				glVertex3f(region.cutting_points[i][0], region.cutting_points[i][1], 0.0);
				glVertex3f(region.cutting_points[i + 1][0], region.cutting_points[i + 1][1], 0.0);
				glEnd();
			}

			glPointSize(point_size);
			glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < region.cutting_points.size(); i++)
			{
				glVertex3f(region.cutting_points[i][0], region.cutting_points[i][1], 0.0);
			}
			glEnd();
		}


		if (draw_save_critical_points)
		{
			glPointSize(point_size);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < region.sdg.save_critical_points.size(); i++)
			{
				glVertex3f(region.sdg.save_critical_points[i][0], region.sdg.save_critical_points[i][1], 0.0);
			}
			glEnd();
		}

		
		if (false)
		{
			for (int j = 0; j < region.sdg.mas.size(); j++)
			{
				if (j != iiiiii)
				{
					continue;
				}

				glPointSize(point_size);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_POINTS);
				if (iiiiii >= 0 && iiiiii<region.sdg.mas.size())
				{
					for (int i = 0; i < region.sdg.mas[iiiiii].size(); i++)
					{
						glVertex3f(region.sdg.medial_axis_points[region.sdg.mas[iiiiii][i]][0], region.sdg.medial_axis_points[region.sdg.mas[iiiiii][i]][1], 0.0);
					}
				}
				glEnd();
			}

		}


		if (draw_minimal_points)
		{
			glPointSize(point_size);
			glColor3f(255.0 / 255.0, 255.0 / 255.0, 26.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < region.sdg.minimal_points.size(); i++)
			{
				glVertex3f(region.sdg.minimal_points[i][0], region.sdg.minimal_points[i][1], 0.0);
			}
			glEnd();
		}

		if (draw_minimal_points)
		{
			glPointSize(point_size);
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < region.sdg.critical_points.size(); i++)
			{
				glVertex3f(region.sdg.critical_points[i][0], region.sdg.critical_points[i][1], 0.0);
			}
			glEnd();
		}

		if (draw_maximal_points)
		{
			glPointSize(point_size);
			glColor3f(0.0, 0.0, 1.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < region.sdg.maximal_points.size(); i++)
			{
				glVertex3f(region.sdg.maximal_points[i][0], region.sdg.maximal_points[i][1], 0.0);
			}
			glEnd();
		}



		if (draw_inner_concave_points)
		{
			glPointSize(point_size);
			glColor3f(0.0, 255.0 / 255.0, 255.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < region.inner_concave_points.size(); i++)
			{
				glVertex3f(region.inner_concave_points[i][0], region.inner_concave_points[i][1], 0.0);
			}
			glEnd();
		}

		if (draw_polygons_entry_exit)
		{
			for (int i = 0; i < region.polygons_entry_exit.size(); i++)
			{
				glPointSize(point_size);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_POINTS);

				glVertex3f(region.polygons_entry_exit[i][0][0], region.polygons_entry_exit[i][0][1], 0.0);
				glVertex3f(region.polygons_entry_exit[i][1][0], region.polygons_entry_exit[i][1][1], 0.0);

				glEnd();
			}
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


			glPointSize(point_size);
			glColor3f(0.0, 1.0, 1.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < turning_points_exit_temp.size(); i++)
			{
				glVertex3f(turning_points_exit_temp[i][0], turning_points_exit_temp[i][1], 0.0);
			}
			glEnd();

			glPointSize(point_size);
			glColor3f(1.0, 1.0, 0.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < turning_points_entry_temp.size(); i++)
			{
				glVertex3f(turning_points_entry_temp[i][0], turning_points_entry_temp[i][1], 0.0);
			}
			glEnd();


		}





		if (false)
		{
			glPointSize(point_size);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_LINE_LOOP);
			if (iiiiii >= 0 && iiiiii < region.polygons.size())
			{
				for (int i = 0; i < region.polygons[iiiiii].size(); i++)
				{
					glVertex3f(region.polygons[iiiiii][i][0], region.polygons[iiiiii][i][1], 0.0);
				}
			}

			glEnd();

		}


	}


}