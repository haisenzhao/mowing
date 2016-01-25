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

namespace hpcg {

	void ToolpathGenerator::mark_domains(CDT& ct,
		CDT::Face_handle start,
		int index,
		std::list<CDT::Edge>& border)
	{
		if (start->info().nesting_level != -1){
			return;
		}
		std::list<CDT::Face_handle> queue;
		queue.push_back(start);

		while (!queue.empty()){
			CDT::Face_handle fh = queue.front();
			queue.pop_front();
			if (fh->info().nesting_level == -1){
				fh->info().nesting_level = index;
				for (int i = 0; i < 3; i++){
					CDT::Edge e(fh, i);
					CDT::Face_handle n = fh->neighbor(i);
					if (n->info().nesting_level == -1){
						if (ct.is_constrained(e)) border.push_back(e);
						else queue.push_back(n);
					}
				}
			}
		}
	}

	void ToolpathGenerator::mark_domains(CDT& cdt)
	{
		for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
			it->info().nesting_level = -1;
		}

		std::list<CDT::Edge> border;
		mark_domains(cdt, cdt.infinite_face(), 0, border);
		while (!border.empty()){
			CDT::Edge e = border.front();
			border.pop_front();
			CDT::Face_handle n = e.first->neighbor(e.second);
			if (n->info().nesting_level == -1){
				mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
			}
		}
	}


	void ToolpathGenerator::Output_Obj_1(std::string path)
	{
		return;
		std::ofstream file(path);

		if (file.is_open())
		{
			//construct two non-intersecting nested polygons  
			Polygon_2 polygon1;
			polygon1.push_back(Point(0, 0));
			polygon1.push_back(Point(2, 0));
			polygon1.push_back(Point(2, 2));
			polygon1.push_back(Point(0, 2));
			Polygon_2 polygon2;
			polygon2.push_back(Point(0.5, 0.5));
			polygon2.push_back(Point(1.5, 0.5));
			polygon2.push_back(Point(1.5, 1.5));
			polygon2.push_back(Point(0.5, 1.5));

			//Insert the polygons into a constrained triangulation
			CDT cdt;
			//cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);
			//cdt.insert_constraint(polygon2.vertices_begin(), polygon2.vertices_end(), true);

			cdt.insert_constraint(contours.outer_boundary().vertices_begin(), contours.outer_boundary().vertices_end(), true);

			for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
			{
				cdt.insert_constraint(hole_iter->vertices_begin(), hole_iter->vertices_end(), true);
			}

			//Mark facets that are inside the domain bounded by the polygon
			mark_domains(cdt);

			double scale = input_toolpath_size / toolpath_size;

			for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
			{
				//faces_id_1.push_back(fit->vertex(2)->index);
				//faces_id_2.push_back(fit->vertex(1)->index);
				//faces_id_3.push_back(fit->vertex(0)->index);
				if (fit->info().in_domain())
				{
					file << "v " << fit->vertex(2)->point().x() * scale << " " << fit->vertex(2)->point().y() * scale << " 0.0" << std::endl;
					file << "v " << fit->vertex(1)->point().x() * scale << " " << fit->vertex(1)->point().y() * scale << " 0.0" << std::endl;
					file << "v " << fit->vertex(0)->point().x() * scale << " " << fit->vertex(0)->point().y() * scale << " 0.0" << std::endl;
				}

			}

			int i = 1;

			for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
			{
				//faces_id_1.push_back(fit->vertex(2)->index);
				//faces_id_2.push_back(fit->vertex(1)->index);
				//faces_id_3.push_back(fit->vertex(0)->index);
				if (fit->info().in_domain())
				{
					file << "f " << i << " " << i + 1 << " " << i + 2 << std::endl;
					i = i + 3;
				}

			}


			int count = 0;
			for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
				fit != cdt.finite_faces_end(); ++fit)
			{
				if (fit->info().in_domain()) ++count;
			}

			std::cout << "There are " << count << " facets in the domain." << std::endl;

		}
		

		/*
		std::ofstream file(path);

		if (file.is_open())
		{

			double scale = input_toolpath_size / toolpath_size;

			CDT cdt;


			Polygon_2 polygon1;
			polygon1.push_back(Point(0, 0));
			polygon1.push_back(Point(2, 0));
			polygon1.push_back(Point(2, 2));
			polygon1.push_back(Point(0, 2));
			Polygon_2 polygon2;

			polygon2.push_back(Point(0.5, 0.5));
			polygon2.push_back(Point(1.5, 0.5));
			polygon2.push_back(Point(1.5, 1.5));
			polygon2.push_back(Point(0.5, 1.5));

			//Insert the polygons into a constrained triangulation
			cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);
			cdt.insert_constraint(polygon2.vertices_begin(), polygon2.vertices_end(), true);


			//cdt.insert_constraint(contours.outer_boundary().vertices_begin(), contours.outer_boundary().vertices_end(), true);

			for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
			{
				//cdt.insert_constraint(hole_iter->vertices_begin(), hole_iter->vertices_end(), true);
			}

			for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
			{
				//faces_id_1.push_back(fit->vertex(2)->index);
				//faces_id_2.push_back(fit->vertex(1)->index);
				//faces_id_3.push_back(fit->vertex(0)->index);
				if (fit->info().in_domain())
				{
					file << "v " << fit->vertex(2)->point().x() * scale << " " << fit->vertex(2)->point().y() * scale << " 0.0" << std::endl;
					file << "v " << fit->vertex(1)->point().x() * scale << " " << fit->vertex(1)->point().y() * scale << " 0.0" << std::endl;
					file << "v " << fit->vertex(0)->point().x() * scale << " " << fit->vertex(0)->point().y() * scale << " 0.0" << std::endl;
				}

			}

			int i = 1;

			for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
			{
				//faces_id_1.push_back(fit->vertex(2)->index);
				//faces_id_2.push_back(fit->vertex(1)->index);
				//faces_id_3.push_back(fit->vertex(0)->index);
				if (fit->info().in_domain())
				{
					file << "f " << i << " " << i + 1 << " " << i + 2 << std::endl;
					i = i + 3;
				}

			}
		}

		*/
	}

	void ToolpathGenerator::Output_Obj(std::string path)
	{
		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		int layer_number = 3;

		std::ofstream file(path);

		double scale = input_toolpath_size / toolpath_size;

		if (file.is_open())
		{
			file << entry_spiral.size() << std::endl;

			for (int i = contour.size() - 1; i >= 0; i--)
			{
				file << "v " << contour[i][0] * scale << " " << contour[i][1] * scale << " 0.0" << std::endl;
			}

			for (int i = 0; i < contour.size(); i++)
			{
				file << "v " << contour[i][0] * scale << " " << contour[i][1] * scale << " " << layer_number*toolpath_size*scale << std::endl;
			}

			
			file << "f ";
			for (int i = 0; i < contour.size(); i++)
			{
				file << i + 1 << " ";
			}
			file << std::endl;

			file << "f ";
			for (int i = 0; i < contour.size(); i++)
			{
				file << i + 1 + contour.size() << " ";
			}


			file << std::endl;

			for (int i = 0; i < contour.size(); i++)
			{
				file << "f " << contour.size() - 1 - i + 1 << " " << contour.size() - 1 - (i + 1) % contour.size() + 1 << " " << (i + 1) % contour.size() + 1 + contour.size() << " " << i + 1 + contour.size() << std::endl;
			}

		}

		file.clear();
		file.close();

		std::vector<Vector2d>().swap(contour);
	}

	void ToolpathGenerator::OutputPathTwoCircles()
	{
		std::ofstream file("D:\\task2\\SLAM\\3DprintingFramework\\GenerateTestData\\square.txt");

		int sub_int_0 = 8;
		int sub_int_1 = 50;

		file << "2" << std::endl;
		file << sub_int_0 << std::endl;
		for (int i = 0; i < sub_int_0; i++)
		{
			file << 10.0*sin(i * 2 * PI / sub_int_0 + PI / 4.0 + PI / 8.0) << " " << 10.0*cos(i * 2 * PI / sub_int_0 + PI / 4.0 + PI / 8.0) << std::endl;
		}

		//file << "6" << std::endl;
		for (int i = 5; i >= 0; i--)
		{
			//file << exit_d_0*sin(i * 2 * PI / 60.0) << " " << exit_d_0*cos(i * 2 * PI / 60.0) << std::endl;
		}
		file << "end" << std::endl;
	}

	void ToolpathGenerator::OutputPath_Hollow(std::string path)
	{
		if (debug_int_0 > 0)
		{
			for (int i = 0; i < smooth_number; i++)
			{
				PolygonSmoothing();
				DirectlyPolygonSmoothing();
			}

			std::vector<Vector2d> contour;
			for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
			{
				contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
			}
			std::vector<Vector2d> offset;

			Circuit::GenerationOffsetWithClipper(contour, debug_int_0*toolpath_size, offset);

			std::ofstream file(path+".txt");

			if (file.is_open())
			{
				file << "2" << std::endl;

				file << contour.size() << std::endl;
				for (int i = 0; i < contour.size(); i++)
				{
					file << -contour[i][0] << " " << -contour[i][1] << std::endl;
				}

				file << offset.size() << std::endl;
				for (int i = offset.size() - 1; i >= 0; i--)
				{
					file << -offset[i][0] << " " << -offset[i][1] << std::endl;
				}
			}
			file.clear();
			file.close();


			Polygon_2 one_contour;

			for (int i = offset.size() - 1; i >= 0; i--)
			{
				one_contour.push_back(Point_2(offset[i][0], offset[i][1]));
			}
			contours.add_hole(one_contour);
		

			//offsets
			std::vector<Vector2d>().swap(offset);

			for (int i = 0; i < debug_int_0; i++)
			{
				Circuit::GenerationOffsetWithClipper(contour, toolpath_size /2.0+ i*toolpath_size, offset);
				offsets.push_back(offset);
				std::vector<Vector2d>().swap(offset);
			}

			Vector2d input_entry_point, input_exit_point;
			Vector2d output_entry_point, output_exit_point;
			Circuit::FindOptimalEntryExitPoints(toolpath_size, offsets[0], input_entry_point, input_exit_point,debug_points);

			//generate Fermat spiral
			//FermatsSpiralTrick(offsets, input_entry_point, input_exit_point, output_entry_point, output_exit_point);


			RichardMethod(offsets, input_entry_point, input_exit_point, output_entry_point, output_exit_point);


			entry_spirals.push_back(entry_spiral);
			exit_spirals.push_back(exit_spiral);

			std::reverse(exit_spiral.begin(), exit_spiral.end());

			pathes.push_back(entry_spiral);
			pathes.push_back(exit_spiral);

			for (int i = 0; i < entry_spiral.size(); i++)
			{
				one_single_path.push_back(entry_spiral[i]);
			}

			for (int i = 0; i < exit_spiral.size(); i++)
			{
				one_single_path.push_back(exit_spiral[i]);
			}

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

			OutputPath(one_single_path, path + ".dat");
			Output_Obj_1(path + ".obj");
		}

	}

	void ToolpathGenerator::OutputPathDirectly(std::vector<Vector2d> &vecs, std::string path)
	{
		std::vector<int> remove_int;

		if (vecs.size() > 2)
		{
			for (int i = 0; i < vecs.size() - 1; i++)
			{
				double d = Strip::Distance(vecs[i], vecs[(i + 1) % vecs.size()]);
				if (d < 0.00001)
				{
					remove_int.push_back((i + 1) % vecs.size());
				}
			}

			for (int i = remove_int.size() - 1; i >= 0; i--)
			{
				vecs.erase(vecs.begin() + remove_int[i]);
			}
		}

		std::ofstream file(path);

		if (file.is_open())
		{
			file << "1" << std::endl;
			file << vecs.size() << std::endl;
			for (int i = 0; i < vecs.size(); i++)
			{
				file << vecs[i][0] << " " << vecs[i][1] << std::endl;
			}
		}
		file.clear();
		file.close();
	}


	void ToolpathGenerator::InputOneSinglePath(std::string path)
	{

		double scale = toolpath_size / input_toolpath_size;

		std::ifstream inputfile(path, std::ios::in);

		int input_int;
		inputfile >> input_int;

		for (int i = 0; i < input_int; i++)
		{
			Vector2d v;
			inputfile >> v[0]>>v[1];
			one_single_path.push_back(Vector2d(v[0] * scale, v[1] * scale));
		}
	}

	void ToolpathGenerator::OutputPath(std::vector<Vector2d> &vecs, std::string path)
	{
		/*
		std::vector<int> remove_int;

		if (vecs.size() > 2)
		{
			for (int i = 0; i < vecs.size() - 1; i++)
			{
				double d = Strip::Distance(vecs[i], vecs[(i + 1) % vecs.size()]);
				if (d < 0.00001)
				{
					remove_int.push_back((i + 1) % vecs.size());
				}
			}

			for (int i = remove_int.size() - 1; i >= 0; i--)
			{
				vecs.erase(vecs.begin() + remove_int[i]);
			}
		}
		*/

		std::ofstream file(path);

		double scale = input_toolpath_size / toolpath_size;

		if (file.is_open())
		{
			file << "1" << std::endl;
			file << vecs.size() << std::endl;
			for (int i = 0; i < vecs.size(); i++)
			{
				if (one_single_path_boundary[i])
					file << vecs[i][0] * scale << " " << vecs[i][1] * scale <<" 1"<< std::endl;
					else
					file << vecs[i][0] * scale << " " << vecs[i][1] * scale <<" 0"<< std::endl;
			}
		}
		file.clear();
		file.close();
	}

	void ToolpathGenerator::Output_tree(std::string path)
	{
		std::ofstream file(path);

		file << "Mark Newman on Sat Jul 22 05:32:16 2006" << std::endl;
		file << "graph" << std::endl;
		file << "[" << std::endl;
		file << "  directed 0" << std::endl;

		for (int i = 0; i < offsets.size(); i++)
		{
			file << "node" << std::endl;
			file << "[" << std::endl;
			file << "id " << i << std::endl;
			file << "label " << i << std::endl;
			file << "]" << std::endl;
		}


		for (int i = 0; i < offset_graph.size(); i = i + 2)
		{
			int index_0 = offset_graph[i];
			int index_1 = offset_graph[i + 1];

			file << "edge" << std::endl;
			file << "[" << std::endl;

			file << "source " << index_0 << std::endl;
			file << "target " << index_1 << std::endl;

			file << "]" << std::endl;
		}

		file << "]" << std::endl;

		file.clear();
		file.close();
	}

	void ToolpathGenerator::Output_tree(std::vector<int> &nodes, std::vector<int> &edges, std::string path)
	{
		std::ofstream file(path);

		file << "Mark Newman on Sat Jul 22 05:32:16 2006" << std::endl;
		file << "graph" << std::endl;
		file << "[" << std::endl;
		file << "  directed 0" << std::endl;


		for (int i = 0; i < nodes.size(); i++)
		{
			file << "node" << std::endl;
			file << "[" << std::endl;
			file << "id " << nodes[i] << std::endl;
			file << "label " << nodes[i] << std::endl;

			file << "]" << std::endl;
		}

		for (int i = 0; i < edges.size(); i = i + 2)
		{
			file << "edge" << std::endl;
			file << "[" << std::endl;

			file << "source " << edges[i] << std::endl;
			file << "target " << edges[i + 1] << std::endl;

			file << "]" << std::endl;
		}

		file << "]" << std::endl;

		file.clear();
		file.close();
	}

}