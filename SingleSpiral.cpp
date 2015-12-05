#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"


namespace hpcg {

	void ToolpathGenerator::ArchinedeanSpiral()
	{
		for (int i = 0; i < 3; i++)
			PolygonSmoothing();

		ComputeOffsets();

		for (int i = 0; i < offsets.size(); i++)
		{
			std::cout << i << " / " << offsets.size() << std::endl;
			double entry_d_1 = ComputeNextTurningPoint(entry_d_0, toolpath_size, i);
			
			SelectOnePartOffset(i, entry_d_0, entry_d_1, entry_spiral);

			if (i != offsets.size() - 1)
				entry_d_0 = FindNearestPoint(i + 1, entry_spiral[entry_spiral.size() - 1]);
		}
	
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