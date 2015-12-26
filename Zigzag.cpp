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

	void ToolpathGenerator::OutputPathTwoCircles()
	{
		std::ofstream file("D:\\task2\\SLAM\\3DprintingFramework\\GenerateTestData\\twocirle.txt");

		file << "2" << std::endl;
		file << "4" << std::endl;
		for (int i = 0; i < 60; i++)
		{
			file << 10.0*sin(i * 2 * PI / 60.0 + PI / 4.0) << " " << 10.0*cos(i * 2 * PI / 60.0 + PI / 4.0) << std::endl;
		}
		file << "60" << std::endl;
		for (int i = 59; i >= 0; i--)
		{
			file << exit_d_0*sin(i * 2 * PI / 60.0) << " " << exit_d_0*cos(i * 2 * PI / 60.0) << std::endl;
		}
		file << "end" << std::endl;
	}

	void ToolpathGenerator::OutputPath(std::string path)
	{
		std::ofstream file(path);

		double scale = 0.4 / toolpath_size;

		if (file.is_open())
		{
			file << "2" << std::endl;
			file << entry_spiral.size() << std::endl;

			for (int i = 0; i < entry_spiral.size(); i++)
			{
				file << entry_spiral[i][0] * scale << " " << entry_spiral[i][1] * scale << std::endl;
			}

			file << exit_spiral.size() << std::endl;

			for (int i = 0; i < exit_spiral.size(); i++)
			{
				file << exit_spiral[i][0] * scale << " " << exit_spiral[i][1] * scale << std::endl;
			}
		}

		file.clear();
		file.close();

		double l = Strip::GetTotalLength(entry_spiral) + Strip::GetTotalLength(exit_spiral);

		std::cout << "Total length: " << l << std::endl;

	}
	void ToolpathGenerator::OutputPath(std::vector<Vector2d> &vecs, std::string path)
	{
		std::ofstream file(path);

		double scale = 0.4 / toolpath_size;

		if (file.is_open())
		{
			for (int i = 0; i < vecs.size(); i++)
			{
				file << vecs[i][0] * scale << " " << vecs[i][1] * scale << std::endl;
			}
		}
		file.clear();
		file.close();
	}


	void ToolpathGenerator::GenerateZigzagForCircle()
	{
		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::vector<std::vector<Vector2d>> offset_t;
		Circuit::GenerateOffsetHole(contours, toolpath_size / 2.0, offset_t);

		int i = 0;

		for (double y = contours.bbox().ymin(); y < contours.bbox().ymax(); y = y + toolpath_size)
		{
			i++;

			Line_2 l2(Point_2(0.0, y), Point_2(1.0, y));

			std::vector<Vector2d> t;

			for (int j = 0; j < offset_t.size(); j++)
			{
				for (int i = 0; i < offset_t[j].size(); i++)
				{
					Point_2 p0(offset_t[j][i][0], offset_t[j][i][1]);
					Point_2 p1(offset_t[j][(i + 1) % offset_t[j].size()][0], offset_t[j][(i + 1) % offset_t[j].size()][1]);

					CGAL::Object result = CGAL::intersection(Segment_2(p0, p1), l2);

					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
					{
						t.push_back(Vector2d(ipoint->x(), ipoint->y()));
					}
				}

			}


			if (t.size() == 4)
			{
				std::vector<Vector2d> vecs;
				SortVector(t, vecs);

				if (i % 2 == 1)
				{
					if (vecs[0][0] < vecs[1][0])
					{
						entry_spiral.push_back(vecs[0]);
						entry_spiral.push_back(vecs[1]);
					}
					else
					{
						entry_spiral.push_back(vecs[1]);
						entry_spiral.push_back(vecs[0]);
					}
				}
				else
				{
					if (vecs[0][0] < vecs[1][0])
					{
						entry_spiral.push_back(vecs[1]);
						entry_spiral.push_back(vecs[0]);
					}
					else
					{
						entry_spiral.push_back(vecs[0]);
						entry_spiral.push_back(vecs[1]);
					}
				}

				if (i % 2 == 1)
				{
					if (vecs[2][0] < vecs[3][0])
					{
						exit_spiral.push_back(vecs[2]);
						exit_spiral.push_back(vecs[3]);
					}
					else
					{
						exit_spiral.push_back(vecs[3]);
						exit_spiral.push_back(vecs[2]);
					}
				}
				else
				{
					if (vecs[2][0] < vecs[3][0])
					{
						exit_spiral.push_back(vecs[3]);
						exit_spiral.push_back(vecs[2]);
					}
					else
					{
						exit_spiral.push_back(vecs[2]);
						exit_spiral.push_back(vecs[3]);
					}
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
	}

	void ToolpathGenerator::GenerateZigzag()
	{
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
	}
}

