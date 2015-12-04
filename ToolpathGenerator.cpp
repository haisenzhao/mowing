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

	std::vector<PixelIndex> insertpixels;


	int iiiiii = 0;
	void  ToolpathGenerator::StepDebug()
	{
		//SearchOnePath(false);
		//SearchOneStepforSpiralExit(false);

		iiiiii++;
		//draw_offset = !draw_offset;


	//	PolygonSmoothing();
	}



	void ToolpathGenerator::Init(TMeshProcessor* render)
	{

		image_space.pixel_size = 0.05;
		image_space.toolpath_size = 0.2;
		toolpath_size = 0.2;

		spiral_affect_index = 0;

		spiral_entry_affect_index = 0;
		spiral_exit_affect_index = 0;

		draw_boundary = true;
		draw_pixels = true;
		draw_aixs =false;
		draw_contour = true;
		draw_offset = true;
		draw_spiral = false;

		linewidth = 5;

		spiral_direction = false;

		first_two_step_bool = false;

		entry_d_0 = 0.8;
		exit_d_0 = 0.2;

		//parameter
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
		//parameter
		
		image_space.toolpath_size = toolpath_size;

		m_render = render;
		LoadContour();
		
		//FermatSpiral();
		//GenerateToolpathBasedonImageSpace();
		//OffsetBasedFermatSpiral();

		OffsetBasedFermatSpiral1();

		//ArchinedeanSpiral();

		//Zigzag();

		/*
		iss = CGAL::create_interior_straight_skeleton_2(contours);

		std::cout << "Straight skeleton with " << iss->size_of_vertices()
			<< " vertices, " << iss->size_of_halfedges()
			<< " halfedges and " << iss->size_of_faces()
			<< " faces" << std::endl;

		*/

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
							if (image_space.pixels[i][j].filled_entry_exit)
							{
								glColor3f(0.8, 0.8, 0.8);
								glBegin(GL_POLYGON);
								DrawPixel(i, j);
								glEnd();
							}
							else
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
		}

		//draw contour
		if (draw_contour)
		{
			glLineWidth(linewidth+1);
			glColor3f(65.0/255.0, 128.0/255.0, 187.0/255.0);

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
		if (draw_offset)
		{
			glLineWidth(2);
			glColor3f(0.5, 0.5, 0.0);

			for (int i = 0; i < toolpath.size(); i++)
			{
				glBegin(GL_LINE_LOOP);
				for (int j = 0; j < toolpath[i].size(); j++)
				{
					glVertex3f(toolpath[i][j][0], toolpath[i][j][1], 0.0);
				}
				glEnd();
			}
		}

		//draw boundary
		if (draw_boundary)
		{
			glLineWidth(3);
			glColor3f(0.0, 0.0, 0.0);

			for (int i = 0; i < image_space.pixels.size(); i++)
			{
				for (int j = 0; j < image_space.pixels[i].size(); j++)
				{
					if (image_space.pixels[i][j].boundary)
					{
						glBegin(GL_POLYGON);
						DrawPixel(i, j);
						glEnd();
					}
				}
			}
		}

		//draw spiral toolpath
		if (draw_spiral)
		{
			glLineWidth(1);
			
			for (int i = 0; i < spirals.size(); i++)
			{
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_LINE_STRIP);
				for (int j = 0; j < spirals[i].size(); j++)
				{
					//DrawPixel(spirals[i][j].x, spirals[i][j].y);
					glVertex3f(image_space.pixels[spirals[i][j].x][spirals[i][j].y].center[0], image_space.pixels[spirals[i][j].x][spirals[i][j].y].center[1], 0.0);
				}

				glEnd();
			}

			for (int i = 0; i < spiral.size(); i++)
			{
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_POLYGON);
				DrawPixel(spiral[i].x, spiral[i].y);
				glEnd();
			}


			if (spiral.size() > 1)
			{
				glLineWidth(1);

				if (image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_index >= 0)
				{
					glColor3f(0.0, 1.0, 0.0);
					glBegin(GL_POLYGON);
					int related_boundary_index = image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_index;
					DrawPixel(spiral[related_boundary_index].x, spiral[related_boundary_index].y);
					glEnd();
				}

				if (image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_x >= 0 &&
					image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_y >= 0)
				{
					glColor3f(0.0, 1.0, 0.0);
					glBegin(GL_POLYGON);
					DrawPixel(image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_x, image_space.pixels[spiral[spiral.size() - 1].x][spiral[spiral.size() - 1].y].related_boundary_y);
					glEnd();
				}

				if (spiral_affect_index < spiral.size())
				{
					glColor3f(0.0, 1.0, 1.0);
					glBegin(GL_POLYGON);
					DrawPixel(spiral[spiral_affect_index].x, spiral[spiral_affect_index].y);
					glEnd();
				}
			}
		}

		if (false)
		{
			glLineWidth(1);
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_POLYGON);
			DrawPixel(247, 30);
			glEnd();
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

		
		glLineWidth(1);
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_POLYGON);
		DrawPixel(entry_point_x, entry_point_y);
		glEnd();
		glBegin(GL_POLYGON);
		DrawPixel(exit_point_x, exit_point_y);
		glEnd();
	

		if (true)
		{
			glLineWidth(1);
			for (int i = 0; i < spiral_entry.size(); i++)
			{
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_POLYGON);
				DrawPixel(spiral_entry[i].x, spiral_entry[i].y);
				glEnd();
			}

			for (int i = 0; i < spiral_exit.size(); i++)
			{
				glColor3f(0.0, 1.0, 0.0);
				glBegin(GL_POLYGON);
				DrawPixel(spiral_exit[i].x, spiral_exit[i].y);
				glEnd();
			}

			//spiral_entry_affect_index

			if (spiral_exit.size() > 0 && spiral_exit_affect_index>0)
			{
				glColor3f(0.7, 0.2, 0.4);
				glBegin(GL_POLYGON);
				DrawPixel(spiral_exit[spiral_exit_affect_index - 1].x, spiral_exit[spiral_exit_affect_index - 1].y);
				glEnd();
			}

		
			if (spiral_entry.size() > 0 && spiral_entry_affect_index>0)
			{
				glBegin(GL_POLYGON);
				DrawPixel(spiral_entry[spiral_entry_affect_index - 1].x, spiral_entry[spiral_entry_affect_index - 1].y);
				glEnd();
			}


			if (spiral_entry.size() > 201)
			{
				glBegin(GL_POLYGON);
				DrawPixel(spiral_entry[200].x, spiral_entry[200].y);
				glEnd();
			}

			if (spiral_exit.size() > 406)
			{
				glBegin(GL_POLYGON);
				DrawPixel(spiral_exit[406].x, spiral_exit[406].y);
				glEnd();
			}

		}
	

		if (false)
		{
			glPointSize(20);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_POINTS);
			glVertex3f(toolpath[0][0].x, toolpath[0][0].y, 0.0);
			glEnd();
		}


		if (!draw_offset)
		{
			glLineWidth(linewidth);
			glColor3f(255.0/255.0, 42.0/255.0, 26.0/255.0);
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < entry_spiral.size() - iiiiii; i++)
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

		if (true)
		{
			glLineWidth(2);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < aaaa.size(); i++)
			{
				glVertex3f(aaaa[i][0], aaaa[i][1], 0.0);
			}
			glEnd();

			glLineWidth(2);
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < bbbb.size(); i++)
			{
				glVertex3f(bbbb[i][0], bbbb[i][1], 0.0);
			}
			glEnd();

			glLineWidth(2);
			glColor3f(0.0, 0.0, 1.0);
			glBegin(GL_LINE_STRIP);

			if (aaaa.size() > 0 && bbbb.size())
			{
				glVertex3f(aaaa[aaaa.size() - 1][0], aaaa[aaaa.size() - 1][1], 0.0);
				glVertex3f(bbbb[bbbb.size() - 1][0], bbbb[bbbb.size() - 1][1], 0.0);
			}
		
			glEnd();
		}


		if (aaaa.size() > 0 && bbbb.size())
		{
			glPointSize(6);
			glColor3f(1.0, 0.0, 0.0);
			glBegin(GL_POINTS);
			glVertex3f(aaaa[0][0], aaaa[0][1], 0.0);
			glEnd();

			glPointSize(6);
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_POINTS);
			glVertex3f(bbbb[0][0], bbbb[0][1], 0.0);
			glEnd();
		}

		if (false)
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

		if (false)
		if (draw_offset)
		{

			glPointSize(linewidth*2.0);
			glColor3f(2 / 255.0, 126 / 255.0, 18.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < dddd.size(); i++)
			{
				glVertex3f(dddd[i][0], dddd[i][1], 0.0);
			}
			glEnd();

			glPointSize(linewidth * 2.0);
			glColor3f(255.0 / 255.0, 42.0 / 255.0, 26.0 / 255.0);
			glBegin(GL_POINTS);
			for (int i = 0; i < cccc.size(); i++)
			{
				glVertex3f(cccc[i][0], cccc[i][1], 0.0);
			}
			glEnd();
		}


		/*
		Halfedge_const_handle null_halfedge;
		Vertex_const_handle   null_vertex;
		for (Halfedge_const_iterator i = iss->halfedges_begin(); i != iss->halfedges_end(); ++i)
		{

			if (i->is_inner_bisector())
			{
				glLineWidth(5);
				glColor3f(1.0, 0.0, 0.0);
				glBegin(GL_LINE_STRIP);

				glVertex3f(i->opposite()->vertex()->point().x(), i->opposite()->vertex()->point().y(), 0.0);
				glVertex3f(i->vertex()->point().x(), i->vertex()->point().y(), 0.0);

				glEnd();
			}
			else
			{
				glLineWidth(5);
				glColor3f(0.0, 0.0, 1.0);
				glBegin(GL_LINE_STRIP);

				
				glVertex3f(i->opposite()->vertex()->point().x(), i->opposite()->vertex()->point().y(), 0.0);
				glVertex3f(i->vertex()->point().x(), i->vertex()->point().y(), 0.0);

				glEnd();
			}
		}
		*/
	}

	void ToolpathGenerator::GenerateToolpath()
	{
		double lOffset = toolpath_size/2.0;
		PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset, contours);
		while (offset_polygons.size()>0)
		{
			for (PolygonPtrVector::const_iterator pi = offset_polygons.begin(); pi != offset_polygons.end(); ++pi)
			{
				std::vector<Vector2d> one_path;
				for (Polygon_2::Vertex_const_iterator vi = (**pi).vertices_begin(); vi != (**pi).vertices_end(); ++vi)
				{
					one_path.push_back(Vector2d((*vi).x(), (*vi).y()));
				}

				toolpath.push_back(one_path);
			}

			lOffset = lOffset + toolpath_size;
			offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset, contours);
		}
	}


	void ToolpathGenerator::GenerateToolpathBasedonOffset()
	{
		std::cout << "GenerateToolpathBasedonOffset" << std::endl;
	}


}