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

	int iiiiii = -1;
	void  ToolpathGenerator::StepDebug()
	{
		iiiiii++;
	}

	void ToolpathGenerator::Init(TMeshProcessor* render)
	{
		image_space.pixel_size = 0.05;
		image_space.toolpath_size = 0.2;
		toolpath_size = 0.2;

		draw_pixels = false;
		draw_aixs =false;
		draw_contour = true;
		draw_offsets = false;
		draw_entry_exit_spiral = true;
		draw_turning_points = true;

		draw_spiral = false;
		draw_voronoi = false;
		draw_medial_axis = false;
		draw_minimal_points = false;
		draw_entry_exit_points = true;
		draw_cutting_points = true;

		linewidth = 5;

		entry_d_0 = 0.8;
		exit_d_0 = 0.2;

		std::string path = "D:\\test.txt";
		std::ifstream file(path, std::ios::in);

		if (!file)
		{
			std::cout << "" << std::endl;
			return;
		}
		file >> linewidth >> input_int_2;
		file >> image_space.pixel_size >> toolpath_size;
		file >> entry_d_0 >> exit_d_0;
		file.clear();
		file.close();
		
		image_space.toolpath_size = toolpath_size;

		m_render = render;
		LoadContour();
		
		
		FillingAlgorithm();

		//FermatSpiral();

		//ArchinedeanSpiral();
		//GenerateZigzag();
		//OutputPath(entry_spiral, "D:\\123.dat");
		
	}
	
	void ToolpathGenerator::LoadContour()
	{
		//load
		std::string path = "D:\\task2\\SLAM\\3DprintingFramework\\GenerateTestData\\test.txt";
		std::ifstream file(path, std::ios::in);
		
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
			glLineWidth(linewidth + 1);
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
			glLineWidth(2);
			glColor3f(0.5, 0.5, 0.0);

			for (int i = 0; i < offsets.size(); i++)
			{
				glBegin(GL_LINE_LOOP);
				for (int j = 0; j < offsets[i].size(); j++)
				{
					glVertex3f(offsets[i][j][0], offsets[i][j][1], 0.0);
				}
				glEnd();
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
			glLineWidth(linewidth);
			glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < entry_spiral.size(); i++)
			{
				glVertex3f(entry_spiral[i][0], entry_spiral[i][1], 0.0);
			}
			glEnd();

			glLineWidth(linewidth);
			glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < exit_spiral.size(); i++)
			{
				glVertex3f(exit_spiral[i][0], exit_spiral[i][1], 0.0);
			}
			glEnd();

			if (true)
			{
				glLineWidth(linewidth);
				glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
				glBegin(GL_LINE_STRIP);

				if (entry_spiral.size() > 0 && exit_spiral.size())
				{
					glVertex3f(entry_spiral[entry_spiral.size() - 1][0], entry_spiral[entry_spiral.size() - 1][1], 0.0);
					glVertex3f(exit_spiral[exit_spiral.size() - 1][0], exit_spiral[exit_spiral.size() - 1][1], 0.0);
				}
				glEnd();
			}
		}



		if (draw_entry_exit_points)
		{
			if (entry_spiral.size() > 0 && exit_spiral.size())
			{
				glPointSize(linewidth*1.5);
				glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
				glBegin(GL_POINTS);
				glVertex3f(entry_spiral[0][0], entry_spiral[0][1], 0.0);
				glEnd();

				glPointSize(linewidth * 1.5);
				glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
				glBegin(GL_POINTS);
				glVertex3f(exit_spiral[0][0], exit_spiral[0][1], 0.0);
				glEnd();
			}
		}

		if (draw_voronoi)
		{
			for (int i = 0; i < voronoi_edge_points.size(); i = i + 2)
			{
				glLineWidth(2);
				glColor3f(0.5, 0.5, 0.0);
				glBegin(GL_LINE_STRIP);
				glVertex3f(voronoi_edge_points[i][0], voronoi_edge_points[i][1], 0.0);
				glVertex3f(voronoi_edge_points[(i + 1) % voronoi_edge_points.size()][0], voronoi_edge_points[(i + 1) % voronoi_edge_points.size()][1], 0.0);
				glEnd();
			}

		}

		if (draw_medial_axis)
		{
			for (int i = 0; i < medial_axis_points.size(); i = i + 2)
			{
				glLineWidth(2 * 2);
				glColor3f(0.4, 0.4, 0.0);
				glBegin(GL_LINE_STRIP);
				glVertex3f(medial_axis_points[i][0], medial_axis_points[i][1], 0.0);
				glVertex3f(medial_axis_points[(i + 1) % medial_axis_points.size()][0], medial_axis_points[(i + 1) % medial_axis_points.size()][1], 0.0);
				glEnd();
			}
			if (false)
				for (int i = 0; i < medial_axis_points.size(); i = i + 2)
				{
					glPointSize(2 * 2 + 2);
					glColor3f(1.0, 0.0, 0.0);
					glBegin(GL_POINTS);
					glVertex3f(medial_axis_points[i][0], medial_axis_points[i][1], 0.0);
					glVertex3f(medial_axis_points[(i + 1) % medial_axis_points.size()][0], medial_axis_points[(i + 1) % medial_axis_points.size()][1], 0.0);
					glEnd();
				}
		}

		if (draw_minimal_points)
		{
			glPointSize(linewidth * 2.0);
			glColor3f(255.0 / 255.0, 255.0 / 255.0, 26.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < minimal_points.size(); i++)
			{
				glVertex3f(minimal_points[i][0], minimal_points[i][1], 0.0);
			}
			glEnd();
		}


		if (draw_polygons_entry_exit)
		{
			for (int i = 0; i < polygons_entry_exit.size(); i++)
			{
				glPointSize(linewidth * 2.0);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_POINTS);

				glVertex3f(polygons_entry_exit[i][0][0], polygons_entry_exit[i][0][1], 0.0);
				glVertex3f(polygons_entry_exit[i][1][0], polygons_entry_exit[i][1][1], 0.0);

				glEnd();
			}
		}

		if (draw_cutting_points)
		{
			for (int i = 0; i < cutting_points.size(); i = i + 2)
			{
				glPointSize(linewidth*2.0);
				glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
				glBegin(GL_LINE_STRIP);
				glVertex3f(cutting_points[i][0], cutting_points[i][1], 0.0);
				glVertex3f(cutting_points[i + 1][0], cutting_points[i + 1][1], 0.0);
				glEnd();
			}

			glPointSize(linewidth*2.0 + 4);
			glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < cutting_points.size(); i++)
			{
				glVertex3f(cutting_points[i][0], cutting_points[i][1], 0.0);
			}
			glEnd();
		}

		if (draw_turning_points)
		{
			glPointSize(linewidth*2.0 + 4);
			glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < turning_points_exit.size(); i++)
			{
				glVertex3f(turning_points_exit[i][0], turning_points_exit[i][1], 0.0);
			}
			glEnd();

			glPointSize(linewidth * 2.0);
			glColor3f(255.0 / 255.0, 255.0 / 255.0, 26.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < turning_points_entry.size(); i++)
			{
				glVertex3f(turning_points_entry[i][0], turning_points_entry[i][1], 0.0);
			}
			glEnd();
		}

	}


}