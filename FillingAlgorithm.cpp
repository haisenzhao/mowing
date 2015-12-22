#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"

#include "Circuit.h"

namespace hpcg {

	int my_cmp(int p1, int  p2)
	{
		return p1 > p2;
	}

	void ToolpathGenerator::FillingAlgorithm()
	{
		for (int i = 0; i < smooth_number; i++)
			PolygonSmoothing();

		//create the region
		std::vector<Vector2d> contour;
		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}
		//ComputeOffsets(contour);

		region = Region(contour);
		std::vector<Vector2d>().swap(contour);

		region.DetectInnerConcavePoints();
		region.sdg.GenerateRegionMedialAxisPoints();
		region.sdg.GenerateVoronoiEdgePoints();
		region.sdg.ComputePointsDegree();
		region.sdg.DetectMaximalAndMinimalPoints();
		region.sdg.DecomposeMedialAxis1();
		//region.sdg.DetectCriticalPoints1(region.inner_concave_points, toolpath_size, toolpath_size / 8.0);
		region.sdg.DetectCriticalPoints(region.inner_concave_points, toolpath_size);
		region.ComputeCuttingPoints();
		region.DecomposeSubregions();

		region.GenerateConnectedGraph(toolpath_size);
		region.ComputeTSPlikePath();
		region.ComputeEntryAndExitPoints();

		if (region.polygons.size() == 1)
		{
			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, region.polygons[0]);
			Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, region.polygons[0]);

			FermatsSpiralTrick(region.polygons[0], entry_p_0, exit_p_0);
			//FermatsSpiralSmooth1(region.polygons[0], entry_p_0, exit_p_0);
			//FermatsSpiralSmooth(region.polygons[0], entry_p_0, exit_p_0);
			//FermatSpiral(region.polygons[0], entry_p_0, exit_p_0);
			region.entry_spirals.push_back(entry_spiral);
			region.exit_spirals.push_back(exit_spiral);
		}
		else
		{
			for (int i = 0; i < region.polygons.size(); i++)
			{
				if (i != 2)
					continue;

				std::vector<std::vector<Vector2d>>().swap(offsets);
				std::vector<Vector2d>().swap(entry_spiral);
				std::vector<Vector2d>().swap(exit_spiral);
				FermatSpiral(region.polygons[i], region.polygons_entry_exit[i][0], region.polygons_entry_exit[i][1]);
				//FermatsSpiralTrick(region.polygons[i], region.polygons_entry_exit[i][0], region.polygons_entry_exit[i][1]);
				//FermatsSpiralSmooth(region.polygons[i], region.polygons_entry_exit[i][0], region.polygons_entry_exit[i][1]);
				region.entry_spirals.push_back(entry_spiral);
				region.exit_spirals.push_back(exit_spiral);

			}
		}

		//Generate Fermat spiral for each sub-region
		//GenerateOffsetsForAllPolygons();
	}
	      
	void ToolpathGenerator::GenerateOffsetsForAllPolygons()
	{
		for (int i = 0; i < region.polygons.size(); i++)
		{
			double lOffset = toolpath_size/2.0;

			do
			{
				std::vector<std::vector<Vector2d>> offset;
				
				Circuit::GenerateOffset(region.polygons[i], lOffset, offset);
				lOffset = lOffset + toolpath_size;
				if (offset.size() > 0)
				{
					for (int j = 0; j < offset.size(); j++)
					{
						offsets.push_back(offset[j]);
						std::vector<Vector2d>().swap(offset[j]);
					}

					std::vector<std::vector<Vector2d>>().swap(offset);
				}
				else
				{
					std::vector<std::vector<Vector2d>>().swap(offset);
					break;
				}

			} while (true);

		}
	}
}