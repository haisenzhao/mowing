#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"


namespace hpcg {


	void ToolpathGenerator::OutputPath(std::vector<Vector2d> &vecs, std::string path)
	{
		std::ofstream file(path);

		if (file.is_open())
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				file << vecs[i][0] << " " << vecs[i][1] << std::endl;
			}
		}
		file.clear();
		file.close();

	}

	void ToolpathGenerator::GenerateZigzag()
	{
		for (int i = 0; i < 4; i++)
			PolygonSmoothing();

		std::vector<Vector2d> offset;
		GenerateOffset(false, -1, toolpath_size / 2.0, offset);

		int i = 0;

		for (double y = contours.bbox().ymin(); y < contours.bbox().ymax(); y = y + toolpath_size)
		{
			i++;

			Line_2 l2(Point_2(0.0, y), Point_2(1.0, y));

			std::vector<Vector2d> t;

			for (int i = 0; i < offset.size(); i++)
			{
				Point_2 p0(offset[i][0], offset[i][1]);
				Point_2 p1(offset[(i + 1) % offset.size()][0], offset[(i + 1) % offset.size()][1]);

				CGAL::Object result = CGAL::intersection(Segment_2(p0, p1), l2);

				if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
				{
					t.push_back(Vector2d(ipoint->x(), ipoint->y()));
				}
			}

			if (t.size() == 2)
			{
				if (i % 2 == 1)
				{
					if (t[0][0] < t[1][0])
					{
						entry_spiral.push_back(t[0]);
						entry_spiral.push_back(t[1]);
					}
					else
					{
						entry_spiral.push_back(t[1]);
						entry_spiral.push_back(t[0]);
					}
				}
				else
				{
					if (t[0][0] < t[1][0])
					{
						entry_spiral.push_back(t[1]);
						entry_spiral.push_back(t[0]);
					}
					else
					{
						entry_spiral.push_back(t[0]);
						entry_spiral.push_back(t[1]);
					}
				}
			}
		}

	}
}

