#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"
#include "Strip.h"
#include "Circuit.h"

namespace hpcg {

	
	void ToolpathGenerator::ArchinedeanSpiralTrick(std::vector<Vector2d> &contour)
	{
		for (int i = 0; i <smooth_number; i++)
			PolygonSmoothing();

		std::vector<Vector2d> contour1;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour1.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::cout << "Contour nb: " << contour1.size() << std::endl;

		Circuit::ComputeOffsets(toolpath_size, contour1, offsets);

		for (int i = 0; i < offsets.size() - 1; i++)
		{
			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offsets[i]);
			double entry_d_1 = Circuit::FindNearestPointPar(entry_p_0, offsets[i + 1]);
			Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, offsets[i + 1]);

			std::vector<Vector2d> entry_half;
			Circuit::SelectOnePartOffset(offsets[i], entry_d_0, entry_half);

			double entry_cutting_d = Strip::IntersectPoint(entry_half, Circuit::GetRelatedLine(entry_p_1, offsets[i + 1]));

			double d0 = Circuit::DeltaDEuclideanDistance(entry_d_0, toolpath_size, offsets[i]);
			Vector2d v = Circuit::GetOnePointFromOffset(d0, offsets[i]);
			d0 = Strip::FindNearestPointPar(v, entry_half);
			if (entry_cutting_d>d0)
			{
				entry_cutting_d = d0;
			}

			std::vector<Vector2d> path;
			Strip::SelectOnePart(entry_half, 0.0, entry_cutting_d, path);

			for (int j = 0; j < path.size(); j++)
			{
				entry_spiral.push_back(path[j]);
			}
			entry_spiral.push_back(entry_p_1);


			std::vector<Vector2d>().swap(path);
			std::vector<Vector2d>().swap(entry_half);

			entry_d_0 = entry_d_1;
		}


		double entry_d_1 = ComputeNextTurningPoint(entry_d_0, toolpath_size, offsets.size() - 1);

		Circuit::SelectOnePartOffset(offsets[offsets.size() - 1], entry_d_0, entry_d_1, entry_spiral);

	}

	void ToolpathGenerator::ArchinedeanSpiralSmooth(std::vector<Vector2d> &contour)
	{
		for (int i = 0; i <smooth_number; i++)
			PolygonSmoothing();

		std::vector<Vector2d> contour1;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour1.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::cout << "Contour nb: " << contour1.size() << std::endl;

		Circuit::ComputeOffsets(toolpath_size, contour1, offsets);

		for (int i = 0; i < offsets.size()-1; i++)
		{
			std::vector<Vector2d> current_half;
			std::vector<Vector2d> next_half;

			Circuit::SelectOnePartOffset(offsets[i], entry_d_0, entry_d_0, current_half);
			for (int j = 0; j < current_half.size(); j++)
			{
				double d = Circuit::FindNearestPointPar(current_half[j], offsets[i + 1]);
				Vector2d v = Circuit::GetOnePointFromOffset(d, offsets[i + 1]);
				next_half.push_back(v);
			}

			Strip::StripDeformation(current_half, next_half);
			
			for (int j = 0; j < current_half.size(); j++)
			{
				entry_spiral.push_back(current_half[j]);
			}

			std::vector<Vector2d>().swap(current_half);
			std::vector<Vector2d>().swap(next_half);

			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offsets[i]);

			entry_d_0 = Circuit::FindNearestPointPar(entry_p_0, offsets[i + 1]);

		}
	}

	void ToolpathGenerator::ArchinedeanSpiral(std::vector<Vector2d> &contour)
	{
		for (int i = 0; i <smooth_number; i++)
			PolygonSmoothing();

		std::vector<Vector2d> contour1;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour1.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::cout << "Contour nb: " << contour1.size() << std::endl;

		Circuit::ComputeOffsets(toolpath_size, contour1, offsets);

		for (int i = 0; i < offsets.size(); i++)
		{
			std::cout << i << " / " << offsets.size() << std::endl;
			double entry_d_1 = ComputeNextTurningPoint(entry_d_0, toolpath_size, i);
			
			Circuit::SelectOnePartOffset(offsets[i], entry_d_0, entry_d_1, entry_spiral);
			Vector2d v = Circuit::GetOnePointFromOffset(entry_d_0,offsets[i]);

			if (i != offsets.size() - 1)
				entry_d_0 = Circuit::FindNearestPointPar(v, offsets[i + 1]);
		}

		if (abs(entry_spiral[entry_spiral.size() - 1][0]) < 0.00001&&abs(entry_spiral[entry_spiral.size() - 1][1]) < 0.00001)
		{
			entry_spiral.erase(entry_spiral.begin() + entry_spiral.size()-1);
		}

	}


	void ToolpathGenerator::DirectlyPolygonSmoothing()
	{
		std::vector<Vector2d> contour;


		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));

		}
		std::vector<Vector2d> smooth_contour;

		for (int i = 0; i < contour.size(); i++)
		{
			Vector2d p0(contour[i][0], contour[i][1]);
			Vector2d p1(contour[(i + 1) % contour.size()][0], contour[(i + 1) % contour.size()][1]);
			Vector2d p2(contour[(i + 2) % contour.size()][0], contour[(i + 2) % contour.size()][1]);

			smooth_contour.push_back(Vector2d((p0[0] + p1[0] + p2[0]) / 3.0, (p0[1] + p1[1] + p2[1]) / 3.0));
		}

		Polygon_2 one_contour;
		for (int i = 0; i < smooth_contour.size(); i++)
		{
			one_contour.push_back(Point_2(smooth_contour[i][0], smooth_contour[i][1]));

		}
		contours.clear();
		contours = Polygon_with_holes(one_contour);

		std::vector<Vector2d>().swap(contour);
		std::vector<Vector2d>().swap(smooth_contour);

	}

	void ToolpathGenerator::PolygonSmoothing()
	{
		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::vector<Vector2d> sub_contour;

		for (int i = 0; i < contour.size(); i++)
		{
			Vector2d p0(contour[i][0], contour[i][1]);
			Vector2d p1(contour[(i + 1) % contour.size()][0], contour[(i + 1) % contour.size()][1]);
			Vector2d p2((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0);
			sub_contour.push_back(p0);
			sub_contour.push_back(p2);
		}

		std::vector<Vector2d>().swap(contour);

		Polygon_2 one_contour;

		for (int i = 0; i < sub_contour.size(); i++)
		{
			if (i % 2 == 1)
			{
				Vector2d p0(sub_contour[i][0], sub_contour[i][1]);
				Vector2d p1(sub_contour[(i + 1) % sub_contour.size()][0], sub_contour[(i + 1) % sub_contour.size()][1]);
				Vector2d p2(sub_contour[(i + 2) % sub_contour.size()][0], sub_contour[(i + 2) % sub_contour.size()][1]);

				sub_contour[(i + 1) % sub_contour.size()][0] = (p0[0] + p1[0] + p2[0]) / 3.0;
				sub_contour[(i + 1) % sub_contour.size()][1] = (p0[1] + p1[1] + p2[1]) / 3.0;
			}
		}

		for (int i = 0; i < sub_contour.size(); i++)
		{
			one_contour.push_back(Point_2(sub_contour[i][0], sub_contour[i][1]));

		}
		std::vector<Vector2d>().swap(sub_contour);
		
		contours.clear();

		contours = Polygon_with_holes(one_contour);

	}
}