#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"
#include "Circuit.h"


namespace hpcg {

	Vector2d my_cmp(Vector2d v0, Vector2d v1)
	{
		if (v0[0] > v1[0])
			return v0;
		else
			return v1;
	}

	int  MinimalVector2d(std::vector<Vector2d> &vecs)
	{
		double min_d = MAXDOUBLE;
		int index = -1;
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i][0] < min_d)
			{
				index = i;
				min_d = vecs[i][0];
			}
		}
		return index;
	}

	void SortVector(std::vector<Vector2d> &vecs_0, std::vector<Vector2d> &vecs_1)
	{

		int nb = vecs_0.size();

		for (int i = 0; i < nb; i++)
		{
			int index = MinimalVector2d(vecs_0);

			vecs_1.push_back(vecs_0[index]);

			vecs_0.erase(vecs_0.begin() + index);
		}

	}

	bool compare_double(const Vector2d &d0, const Vector2d &d1)
	{
		return d0[0] > d1[0];
	}

	void ToolpathGenerator::GenerateZigzagForArbitraryShape()
	{
		if (debug_int_0>0)
			OptimalDirection();

		for (int i = 0; i < smooth_number; i++)
			PolygonSmoothing();

		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::vector<Vector2d> offset;
		Circuit::GenerateOffset(contour, toolpath_size / 2.0, offset);

		int iiii = 0;

		std::vector<std::vector<Vector2d>> segmentses;

		for (double y = contours.bbox().ymin(); y < contours.bbox().ymax(); y = y + toolpath_size)
		{
			iiii++;

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

			std::sort(t.begin(), t.end(), compare_double);

			//contours.outer_boundary().bounded_side

			std::vector<Vector2d> one_segments;


			if (t.size()>1)
			{
				for (int i = 0; i < t.size() - 1; i++)
				{
					Vector2d center((t[i][0] + t[(i + 1) % t.size()][0]) / 2.0, (t[i][1] + t[(i + 1) % t.size()][1]) / 2.0);

					if (contours.outer_boundary().bounded_side(Point_2(center[0], center[1])) == CGAL::ON_BOUNDED_SIDE)
					{
						one_segments.push_back(t[i]);
						one_segments.push_back(t[(i + 1) % t.size()]);
					}
				}
			}
			if (one_segments.size()!=0)
			segmentses.push_back(one_segments);
		}



		std::vector<int> edges;
		std::vector<Vector2d> segments;

		for (int i = 0; i < segmentses.size(); i++)
		{
			for (int j = 0;j< segmentses[i].size(); j++)
			{
				segments.push_back(segmentses[i][j]);
				//entry_spiral.push_back(segmentses[i][j]);
			}
		}

		std::ifstream file("D:\\123456.txt", std::ios::in);


		std::vector<std::vector<int>> int_vecs;


		int int0;
		file >> int0;
		
		for (int i = 0; i < int0; i++)
		{
			int int1;
			file >> int1;
			
			std::vector<int> one_int;

			for (int j = 0; j < int1; j++)
			{
				int int2;
				file >> int2;
				one_int.push_back(int2);
			}

			int_vecs.push_back(one_int);
		}


		file >> debug_int_1;
		file >> debug_int_2;


		for (int i = 0; i < int_vecs[0].size(); i++)
		{
			if (i % 2 == 0)
			{
				entry_spiral.push_back(segments[2 * int_vecs[0][i]]);
				entry_spiral.push_back(segments[2 * int_vecs[0][i] + 1]);
			}
			else
			{
				entry_spiral.push_back(segments[2 * int_vecs[0][i] + 1]);
				entry_spiral.push_back(segments[2 * int_vecs[0][i]]);
			}
		}


		std::vector<Vector2d>().swap(contour);

		for (int i = 0; i < entry_spiral.size(); i++)
		{
			one_single_path.push_back(entry_spiral[i]);
		}

	}

	void ToolpathGenerator::GenerateZigzag()
	{
		if (debug_int_0>0)
		OptimalDirection();

		for (int i = 0; i < smooth_number; i++)
			PolygonSmoothing();

		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::vector<Vector2d> offset;
		Circuit::GenerateOffset(contour, toolpath_size / 2.0, offset);

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
		std::vector<Vector2d>().swap(contour);

		for (int i = 0; i < entry_spiral.size(); i++)
		{
			one_single_path.push_back(entry_spiral[i]);
		}
	}
}

