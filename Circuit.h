#pragma once

#include <iostream>
#include <fstream>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>

#include <ToolpathGenerator.h>
#include <Strip.h>

#include "clipper/clipper.hpp"

namespace hpcg {

	class Circuit
	{
	public:
		Circuit()
		{
		}
		~Circuit()
		{
		}

		//basic functions
		static double Angle2PI(Vector2d a, Vector2d b)
		{
			double ab, a1, b1, cosr;
			ab = a[0] * b[0] + a[1] * b[1];
			a1 = sqrt(a[0] * a[0] + a[1] * a[1]);
			b1 = sqrt(b[0] * b[0] + b[1] * b[1]);
			cosr = ab / a1 / b1;
			double angle = acos(cosr);

			if (a[0] * b[1] - a[1] * b[0] > 0.0)
			{
				angle = 2 * PI - angle;
			}
			return  angle;
		}
		static double AnglePI(Vector2d a, Vector2d b)
		{
			double ab, a1, b1, cosr;
			ab = a[0] * b[0] + a[1] * b[1];
			a1 = sqrt(a[0] * a[0] + a[1] * a[1]);
			b1 = sqrt(b[0] * b[0] + b[1] * b[1]);
			cosr = ab / a1 / b1;
			double angle = acos(cosr);
			return  angle;
		}

		//multiplication cross
		static bool MultiCross(Vector2d a, Vector2d b)
		{
			return (a[0] * b[1] - a[1] * b[0]) > 0.0;
		}

		//check position relation
		static bool CheckEnclosed(std::vector<Vector2d> &input_points_0, std::vector<Vector2d> &input_points_1)
		{
			bool b0 = CheckInside(input_points_0[0], input_points_1);
			bool b1 = CheckInside(input_points_1[0], input_points_0);

			return b0 || b1;
		}

		static bool CheckInside(Vector2d v, std::vector<Vector2d> &input_points)
		{
			Polygon_2 poly;

			for (int i = 0; i < input_points.size(); i++)
			{
				poly.push_back(Point_2(input_points[i][0], input_points[i][1]));
			}

			return poly.bounded_side(Point_2(v[0], v[1])) == CGAL::ON_BOUNDED_SIDE;
		}

		static double GetTotalLength(std::vector<Vector2d> &input_points)
		{
			double length = 0.0;
			for (int i = 0; i < input_points.size(); i++)
			{
				length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
			}
			return length;
		}

		//distance
		static double Distance(Vector2d v, std::vector<Vector2d> &input_points)
		{
			if (input_points.size() <= 1)
			{
				if (input_points.size() == 0)
				{
					return MAXDOUBLE;

				}
				else
				{
					return Strip::Distance(v, input_points[0]);
				}
			}
			else
			{
				double m_d = MAXDOUBLE;
				for (int i = 0; i < input_points.size(); i++)
				{
					double d = sqrt((double)CGAL::squared_distance(Point_2(v.x, v.y), Segment_2(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y))));
					m_d = min(m_d, d);
				}
				return m_d;
			}
		}

		static double Distance(Vector2d v, std::vector<std::vector<Vector2d> > &input_pointsese)
		{
			double min_d = MAXDOUBLE;
			for (int i = 0; i < input_pointsese.size(); i++)
			{
				double d = Distance(v, input_pointsese[i]);
				min_d = std::min(min_d, d);
			}
			return min_d;
		}

		static double Distance(Vector2d v0, Vector2d v1, std::vector<Vector2d> &input_points)
		{
			double m_d = MAXDOUBLE;

			for (int i = 0; i < input_points.size(); i++)
			{
				//double d = sqrt((double)CGAL::squared_distance(Line_2(Point_2(v0.x, v0.y), Point_2(v1.x, v1.y)), Line_2(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y))));
				double d = sqrt((double)CGAL::squared_distance(Segment_2(Point_2(v0.x, v0.y), Point_2(v1.x, v1.y)), Segment_2(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y))));
				m_d = min(m_d, d);
			}
			return m_d;

		}

		static bool DistanceDouble(std::vector<Vector2d> &input_points_0, std::vector<Vector2d> &input_points_1, double toolpath_size)
		{
			bool b = false;

			for (int i = 0; i < input_points_0.size() && !b; i++)
			{
				double d = Distance(input_points_0[i], input_points_0[(i + 1) % input_points_0.size()], input_points_1);

				if (abs(d - toolpath_size) < 0.001)
				{
					b = true;
					return true;
				}
			}
			for (int i = 0; i < input_points_1.size() && !b; i++)
			{
				double d = Distance(input_points_1[i], input_points_1[(i + 1) % input_points_1.size()], input_points_0);

				if (abs(d - toolpath_size) < 0.001)
				{
					b = true;
					return true;
				}
			}

			return b;
		}

		static double Distance(std::vector<Vector2d> &input_points_0, std::vector<Vector2d> &input_points_1)
		{
			double m_d = MAXDOUBLE;

			for (int i = 0; i < input_points_0.size(); i++)
			{
				double d = Distance(input_points_0[i], input_points_0[(i + 1) % input_points_0.size()], input_points_1);
				m_d = min(m_d, d);
			}
			for (int i = 0; i < input_points_1.size(); i++)
			{
				double d = Distance(input_points_1[i], input_points_1[(i + 1) % input_points_1.size()], input_points_0);
				m_d = min(m_d, d);
			}
			return m_d;
		}


		//offsets
		static void GenerateOffsetForOutside(std::vector<Vector2d> &input_points, double d, std::vector<Vector2d> &offset)
		{
			Polygon_2 polygon;

			for (int i = 0; i < input_points.size(); i++)
			{
				polygon.push_back(Point_2(input_points[i][0], input_points[i][1]));
			}

			if (polygon.is_simple())
			{
				if (!polygon.is_clockwise_oriented())
				{
					polygon.reverse_orientation();
				}

				PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(d, polygon);

				for (PolygonPtrVector::const_iterator pi = offset_polygons.begin(); pi != offset_polygons.end(); ++pi)
				{
					for (Polygon_2::Vertex_const_iterator vi = (**pi).vertices_begin(); vi != (**pi).vertices_end(); ++vi)
					{
						offset.push_back(Vector2d((*vi).x(), (*vi).y()));
					}
				}
			}
			else
			{
				assert(false);
			}

		}

		static void GenerateOffset(std::vector<Vector2d> &input_points, double d, std::vector<Vector2d> &offset)
		{
			Polygon_2 polygon;

			for (int i = 0; i < input_points.size(); i++)
			{
				polygon.push_back(Point_2(input_points[i][0], input_points[i][1]));
			}

			if (polygon.is_simple())
			{
				if (polygon.is_clockwise_oriented())
				{
					polygon.reverse_orientation();
				}

				PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(d, polygon);

				for (PolygonPtrVector::const_iterator pi = offset_polygons.begin(); pi != offset_polygons.end(); ++pi)
				{
					for (Polygon_2::Vertex_const_iterator vi = (**pi).vertices_begin(); vi != (**pi).vertices_end(); ++vi)
					{
						offset.push_back(Vector2d((*vi).x(), (*vi).y()));
					}
				}
			}
			else
			{
				assert(false);
			}

		}

		static void GenerateOffset(std::vector<Vector2d> &input_points, double d, std::vector<std::vector<Vector2d>> &offsets)
		{
			Polygon_2 polygon;

			for (int i = 0; i < input_points.size(); i++)
			{
				polygon.push_back(Point_2(input_points[i][0], input_points[i][1]));
			}

			if (polygon.is_simple())
			{
				if (polygon.is_clockwise_oriented())
				{
					polygon.reverse_orientation();
				}

				PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(d, polygon);

				for (PolygonPtrVector::const_iterator pi = offset_polygons.begin(); pi != offset_polygons.end(); ++pi)
				{
					std::vector<Vector2d> offset;

					for (Polygon_2::Vertex_const_iterator vi = (**pi).vertices_begin(); vi != (**pi).vertices_end(); ++vi)
					{
						offset.push_back(Vector2d((*vi).x(), (*vi).y()));
					}
					offsets.push_back(offset);
					std::vector<Vector2d>().swap(offset);
				}
			}
			else
			{
				assert(false);
			}
		}

		static void GenerateOffsetWithClipper(std::vector<Vector2d> &input_points, double d, std::vector<std::vector<Vector2d>> &offsets)
		{

			double scale = 1000000.0;

			ClipperLib::Path subj;
			ClipperLib::Paths solution;

			for (int i = 0; i < input_points.size(); i++)
			{
				subj << ClipperLib::IntPoint(input_points[i][0] * scale, input_points[i][1] * scale);
			}

			ClipperLib::ClipperOffset co;
			co.ArcTolerance = co.ArcTolerance*scale / 1000.0;
			co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
			co.Execute(solution, -d*scale);

			for (int i = 0; i < solution.size(); i++)
			{
				std::vector<Vector2d> offset;

				for (int j = 0; j < solution[i].size(); j++)
				{
					offset.push_back(Vector2d(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
				}
				RemoveClosePoints(offset);

				offsets.push_back(offset);

				std::vector<Vector2d>().swap(offset);
			}
		}

		static void GenerateOffsetForOutside(std::vector<Vector2d> &input_points, double d, std::vector<std::vector<Vector2d>> &offsets)
		{
			Polygon_2 polygon;

			for (int i = 0; i < input_points.size(); i++)
			{
				polygon.push_back(Point_2(input_points[i][0], input_points[i][1]));
			}

			if (polygon.is_simple())
			{
				if (!polygon.is_clockwise_oriented())
				{
					polygon.reverse_orientation();
				}

				PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(d, polygon);

				for (PolygonPtrVector::const_iterator pi = offset_polygons.begin(); pi != offset_polygons.end(); ++pi)
				{
					std::vector<Vector2d> offset;

					for (Polygon_2::Vertex_const_iterator vi = (**pi).vertices_begin(); vi != (**pi).vertices_end(); ++vi)
					{
						offset.push_back(Vector2d((*vi).x(), (*vi).y()));
					}
					offsets.push_back(offset);
					std::vector<Vector2d>().swap(offset);
				}
			}
			else
			{
				assert(false);
			}
		}

		static void GenerateOffsetForOutsideWithClipper(std::vector<Vector2d> &input_points, double d, std::vector<std::vector<Vector2d>> &offsets)
		{
			//using namespace ClipperLib;

			double scale = 1000000.0;

			ClipperLib::Path subj;
			ClipperLib::Paths solution;

			for (int i = 0; i < input_points.size(); i++)
			{
				subj << ClipperLib::IntPoint(input_points[i][0] * scale, input_points[i][1] * scale);
			}

			ClipperLib::ClipperOffset co;
			co.ArcTolerance = co.ArcTolerance*scale / 1000.0;
			co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
			co.Execute(solution, d*scale);

			for (int i = 0; i < solution.size(); i++)
			{
				std::vector<Vector2d> offset;

				for (int j = 0; j < solution[i].size(); j++)
				{
					offset.push_back(Vector2d(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
				}
				RemoveClosePoints(offset);

				offsets.push_back(offset);

				std::vector<Vector2d>().swap(offset);
			}
		}

		static void GenerateOffsetHole(Polygon_with_holes polygon_with_hole, double d, std::vector<std::vector<Vector2d>> &offsets)
		{
			PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(d, polygon_with_hole);

			for (PolygonPtrVector::const_iterator pi = offset_polygons.begin(); pi != offset_polygons.end(); ++pi)
			{
				std::vector<Vector2d> offset;

				Polygon_2 poly_2;

				for (Polygon_2::Vertex_const_iterator vi = (**pi).vertices_begin(); vi != (**pi).vertices_end(); ++vi)
				{
					offset.push_back(Vector2d((*vi).x(), (*vi).y()));
					poly_2.push_back(Point_2((*vi).x(), (*vi).y()));
				}

				if (poly_2.is_clockwise_oriented())
				{
					std::reverse(offset.begin(), offset.end());
				}


				offsets.push_back(offset);
				std::vector<Vector2d>().swap(offset);
			}
		}

		static void GenerationOffsetWithClipper(std::vector<Vector2d> &input_points, double d, std::vector<Vector2d> &offset)
		{
			//using namespace ClipperLib;

			double scale = 1000000.0;

			ClipperLib::Path subj;
			ClipperLib::Paths solution;

			for (int i = 0; i < input_points.size(); i++)
			{
				subj << ClipperLib::IntPoint(input_points[i][0] * scale, input_points[i][1] * scale);
			}

			ClipperLib::ClipperOffset co;
			co.ArcTolerance = co.ArcTolerance*scale / 1000.0;
			co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
			co.Execute(solution, -d*scale);

			for (int i = 0; i < solution.size(); i++)
			{
				for (int j = 0; j < solution[i].size(); j++)
				{
					offset.push_back(Vector2d(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
				}
				RemoveClosePoints(offset);
				break;
			}
		}


		static void RemoveClosePoints(std::vector<Vector2d> &input_points)
		{
			std::vector<int> remove_int;

			if (input_points.size() > 2)
			{
				for (int i = 0; i < input_points.size() - 1; i++)
				{
					double d = Strip::Distance(input_points[i], input_points[(i + 1) % input_points.size()]);
					if (d < 0.00001)
					{
						remove_int.push_back((i + 1) % input_points.size());
					}
				}

				for (int i = remove_int.size() - 1; i >= 0; i--)
				{
					input_points.erase(input_points.begin() + remove_int[i]);
				}
			}
		}

		static void UniformDirection(std::vector<Vector2d> &input_points)
		{
			Polygon_2 polygon_2;

			for (int i = 0; i < input_points.size(); i++)
			{
				polygon_2.push_back(Point_2(input_points[i][0], input_points[i][1]));
			}

			if (polygon_2.is_clockwise_oriented())
			{
				std::reverse(input_points.begin(), input_points.end());
			}
		}

		static void GenerationOffsetsHoleWithClipper(Polygon_with_holes polygon_with_hole, double d, std::vector<std::vector<Vector2d>> &offsets)
		{
			//using namespace ClipperLib;

			double scale = 1000000.0;

			ClipperLib::ClipperOffset co;
			co.ArcTolerance = co.ArcTolerance*scale / 1000.0;

			ClipperLib::Path subj;
			ClipperLib::Paths solution;
			//subj << ClipperLib::IntPoint(0, 1) << ClipperLib::IntPoint(1, 1) << ClipperLib::IntPoint(1, 0) << ClipperLib::IntPoint(0, 0);
			for (Polygon_2::Vertex_iterator ver_iter = polygon_with_hole.outer_boundary().vertices_begin(); ver_iter != polygon_with_hole.outer_boundary().vertices_end(); ver_iter++)
			{
				//contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				subj << ClipperLib::IntPoint(ver_iter->x()*scale, ver_iter->y()*scale);
			}

			co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);

			for (Polygon_with_holes::Hole_iterator hole_iter = polygon_with_hole.holes_begin(); hole_iter != polygon_with_hole.holes_end(); hole_iter++)
			{
				subj.clear();
				for (Polygon_2::Vertex_iterator ver_iter = hole_iter->vertices_begin(); ver_iter != hole_iter->vertices_end(); ver_iter++)
				{
					subj << ClipperLib::IntPoint(ver_iter->x()*scale, ver_iter->y()*scale);
				}
				co.AddPath(subj, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
			}

			co.Execute(solution, -d*scale);

			for (int i = 0; i < solution.size(); i++)
			{
				std::vector<Vector2d> one_path;

				Polygon_2 poly_2;
				for (int j = 0; j < solution[i].size(); j++)
				{
					one_path.push_back(Vector2d(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));

					poly_2.push_back(Point_2(((double)solution[i][j].X) / scale, ((double)solution[i][j].Y) / scale));
				}

				if (poly_2.is_clockwise_oriented())
				{
					std::reverse(one_path.begin(), one_path.end());
				}


				//remove very close points
				RemoveClosePoints(one_path);

				offsets.push_back(one_path);
			}
		}

		static double MaxLengthAxis(std::vector<Vector2d> &input_points)
		{
			double min_x = MAXDOUBLE;
			double min_y = MAXDOUBLE;
			double max_x = -MAXDOUBLE;
			double max_y = -MAXDOUBLE;

			for (int i = 0; i < input_points.size(); i++)
			{
				min_x = min(min_x, (double)input_points[i][0]);
				min_y = min(min_y, (double)input_points[i][1]);
				max_x = max(max_x, (double)input_points[i][0]);
				max_y = max(max_y, (double)input_points[i][1]);
			}

			if (max_x - min_x > max_y - min_y)
			{
				return max_x - min_x;
			}
			else
			{
				return max_y - min_y;
			}
		}
		static double MinLengthAxis(std::vector<Vector2d> &input_points)
		{
			double min_x = MAXDOUBLE;
			double min_y = MAXDOUBLE;
			double max_x = -MAXDOUBLE;
			double max_y = -MAXDOUBLE;

			for (int i = 0; i < input_points.size(); i++)
			{
				min_x = min(min_x, (double)input_points[i][0]);
				min_y = min(min_y, (double)input_points[i][1]);
				max_x = max(max_x, (double)input_points[i][0]);
				max_y = max(max_y, (double)input_points[i][1]);
			}

			if (max_x - min_x > max_y - min_y)
			{
				return max_y - min_y;
			}
			else
			{
				return max_x - min_x;
			}
		}

		static void RemoveMinimalOffsets(std::vector<std::vector<Vector2d>> &offsets, double toolpath_size)
		{
			std::vector<int> remove_int;

			for (int i = 0; i < offsets.size(); i++)
			{
				if (offsets[i].size() <= 2)
				{
					remove_int.push_back(i);
				}
				else
				{
					if (MinLengthAxis(offsets[i])<toolpath_size)
					{
						std::vector<Vector2d> offset;
						GenerationOffsetWithClipper(offsets[i], toolpath_size / 4.0, offset);
						if (offset.size() == 0)
						{
							remove_int.push_back(i);
						}
					}
				}
			}

			for (int i = remove_int.size() - 1; i >= 0; i--)
			{
				offsets.erase(offsets.begin() + remove_int[i]);
			}
		}

		static void ComputeOffsets(double toolpath_size, Polygon_with_holes &contours,
			std::vector<std::vector<Vector2d>> &offsets, std::vector<std::vector<std::vector<Vector2d>>> &offsetses)
		{
			double lOffset = toolpath_size / 2.0;
			std::vector<std::vector<Vector2d>> one_pathes;
			//Circuit::GenerateOffsetHole(contours, lOffset, one_pathes);

			GenerationOffsetsHoleWithClipper(contours, lOffset, one_pathes);

			RemoveMinimalOffsets(one_pathes, toolpath_size);
			while (one_pathes.size() > 0)
			{
				std::cout << "Offsets index: " << offsets.size() << std::endl;

				for (int i = 0; i < one_pathes.size(); i++)
				{
					offsets.push_back(one_pathes[i]);
				}

				offsetses.push_back(one_pathes);

				for (int i = 0; i < one_pathes.size(); i++)
				{
					std::vector<Vector2d>().swap(one_pathes[i]);
				}
				std::vector<std::vector<Vector2d>>().swap(one_pathes);

				lOffset = lOffset + toolpath_size;
				GenerationOffsetsHoleWithClipper(contours, lOffset, one_pathes);
				RemoveMinimalOffsets(one_pathes, toolpath_size);
			}

			//add small offsets
			lOffset = lOffset - toolpath_size / 2.0;
			std::vector<std::vector<Vector2d>>().swap(one_pathes);

			GenerationOffsetsHoleWithClipper(contours, lOffset, one_pathes);
			RemoveMinimalOffsets(one_pathes, toolpath_size);

		}

		static void ComputeOffsets(double toolpath_size, std::vector<Vector2d> &contour, std::vector<std::vector<Vector2d>> &offsets)
		{
			double lOffset = toolpath_size / 2.0;
			std::vector<std::vector<Vector2d>> one_pathes;
			Circuit::GenerateOffset(contour, lOffset, one_pathes);

			while (one_pathes.size() > 0)
			{
				std::cout << "Offsets index: " << offsets.size() << std::endl;

				for (int i = 0; i < one_pathes.size(); i++)
				{
					offsets.push_back(one_pathes[i]);
				}

				for (int i = 0; i < one_pathes.size(); i++)
				{
					std::vector<Vector2d>().swap(one_pathes[i]);
				}
				std::vector<std::vector<Vector2d>>().swap(one_pathes);

				lOffset = lOffset + toolpath_size;
				Circuit::GenerateOffset(contour, lOffset, one_pathes);
			}
		}

		//par system
		static double FindNearestPointPar(Vector2d v, std::vector<Vector2d> &input_points)
		{
			Vector2d n_p;

			double total_length = 0.0;
			for (int i = 0; i < input_points.size(); i++)
			{
				total_length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
			}

			double min_d = MAXDOUBLE;
			int min_i = -1;
			for (int i = 0; i < input_points.size(); i++)
			{
				Point_2 p0(input_points[i].x, input_points[i].y);
				Point_2 p1(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y);

				double l = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), Segment_2(p0, p1)));

				if (l < min_d)
				{
					min_d = l;
					min_i = i;
				}
			}

			if (min_i >= 0)
			{
				double length = 0.0;
				for (int i = 0; i < min_i; i++)
				{
					length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
				}

				Point_2 p0(input_points[min_i].x, input_points[min_i].y);
				Point_2 p1(input_points[(min_i + 1) % input_points.size()].x, input_points[(min_i + 1) % input_points.size()].y);

				double l0 = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), p0));
				double l1 = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), p1));

				if (min_d<l0 &&min_d<l1&&abs(min_d - l0)>0.0001&&abs(min_d - l1)>0.0001)
				{
					double l = sqrt((double)CGAL::squared_distance(p0, p1));
					if (l < 0.00001)
					{
						v[0] = p0[0];
						v[1] = p0[1];
					}
					else
					{
						Vector2d vec(p1[0] - p0[0], p1[1] - p0[1]);
						Vector2d r_vec(vec[1], -vec[0]);
						if (vec[0] < 0.00001)
						{
							r_vec[0] = -vec[1];
							r_vec[1] = vec[0];
						}

						if (vec[1] < 0.00001)
						{
							r_vec[0] = vec[1];
							r_vec[1] = -vec[0];
						}

						CGAL::Object result = CGAL::intersection(Segment_2(p0, p1), Line_2(Point_2(v[0], v[1]), Point_2(v[0] + r_vec[0], v[1] + r_vec[1])));

						if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
						{
							n_p[0] = ipoint->x();
							n_p[1] = ipoint->y();
						}
						else
						{
							assert(false);
						}
					}
				}
				else
				{
					if (l0 < l1)
					{
						n_p[0] = p0[0];
						n_p[1] = p0[1];
					}
					else
					{
						n_p[0] = p1[0];
						n_p[1] = p1[1];
					}

				}

				length += sqrt((double)CGAL::squared_distance(Point_2(input_points[min_i].x, input_points[min_i].y), Point_2(n_p[0], n_p[1])));

				return length / total_length;
			}
			else
			{
				assert(false);
			}

			return -1.0;



		}

		static Vector2d GetOnePointFromOffset(double d, std::vector<Vector2d> &input_points)
		{
			Vector2d v;
			double length = 0.0;
			for (int i = 0; i < input_points.size(); i++)
			{
				length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
			}

			double total_length = length;
			length = 0.0;

			for (int i = 0; i < input_points.size(); i++)
			{
				double l = sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));

				if (d*total_length >= length&&d*total_length <= length + l)
				{
					double ll = (d - length / total_length)*total_length / l;
					v[0] = input_points[i].x + (input_points[(i + 1) % input_points.size()].x - input_points[i].x)*ll;
					v[1] = input_points[i].y + (input_points[(i + 1) % input_points.size()].y - input_points[i].y)*ll;
					break;
				}
				length += l;
			}

			return v;
		}

		static void SelectOnePartOffset(std::vector<Vector2d> &input_points, double d, std::vector<Vector2d> &vecs)
		{
			double d1 = d + 0.00001;
			if (d1 > 1.0)
				d1 = d1 - 1.0;

			SelectOnePartOffset(input_points, d, d1, vecs);
			vecs.erase(vecs.begin() + vecs.size() - 1);
			vecs.push_back(vecs[0]);
		}

		static void SelectOnePartOffset(std::vector<Vector2d> &input_points, double d0, double d1, std::vector<Vector2d> &vecs)
		{
			if (abs(d0 - d1) < 0.0000001)
			{
				Vector2d v = GetOnePointFromOffset(d0, input_points);
				vecs.push_back(v);
			}
			else
			{
				Vector2d v = GetOnePointFromOffset(d0, input_points);
				vecs.push_back(v);

				if (d1 >= 0 && d1 <= 1.0)
				{
					double total_length = GetTotalLength(input_points);

					std::vector<double> vec_ds;

					vec_ds.push_back(0.0);
					double length = 0.0;
					for (int i = 0; i < input_points.size(); i++)
					{
						length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
						vec_ds.push_back(length / total_length);
					}


					if (d0 > d1)
					{
						for (int i = vec_ds.size() - 1; i >= 0; i--)
						{
							if (vec_ds[i]<d0&&vec_ds[i]>d1)
							{
								v = GetOnePointFromOffset(vec_ds[i], input_points);

								if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
								{
									vecs.push_back(v);
								}
							}

							if (vec_ds[i] < d1)
							{
								break;
							}
						}
					}
					else
					{
						for (int i = vec_ds.size() - 1; i >0; i--)
						{
							if (vec_ds[i] < d0)
							{
								v = GetOnePointFromOffset(vec_ds[i], input_points);

								if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
								{
									vecs.push_back(v);
								}
							}
						}

						for (int i = vec_ds.size() - 1; i >0; i--)
						{
							if (vec_ds[i] > d1)
							{
								v = GetOnePointFromOffset(vec_ds[i], input_points);
								if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
								{
									vecs.push_back(v);
								}
							}
						}
					}

					v = GetOnePointFromOffset(d1, input_points);
					if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
					{
						vecs.push_back(v);
					}

					if (abs(vecs[1][0] - vecs[0][0]) < 0.000001&&abs(vecs[1][1] - vecs[0][1]) < 0.000001)
					{
						vecs.erase(vecs.begin());
					}

					std::vector<double>().swap(vec_ds);
				}
			}
		}

		//points moving on circuit
		static double DeltaDEuclideanDistance(double d, double distance, std::vector<Vector2d> &input_points, std::vector<Vector2d> &one_path)
		{
			Vector2d v = Circuit::GetOnePointFromOffset(d, input_points);

			std::vector<Vector2d> vecs;

			int divided_nb = 50;
			for (int i = 0; i < divided_nb; i++)
			{
				Point_2 p0(v[0] + abs(distance)*sin(i * 2 * PI / (double)divided_nb), v[1] + abs(distance)*cos(i * 2 * PI / (double)divided_nb));
				Point_2 p1(v[0] + abs(distance)*sin((i + 1) * 2 * PI / (double)divided_nb), v[1] + abs(distance)*cos((i + 1) * 2 * PI / (double)divided_nb));

				one_path.push_back(Vector2d(p0[0], p0[1]));

				for (int j = 0; j < input_points.size(); j++)
				{
					Point_2 p2(input_points[j].x, input_points[j].y);
					Point_2 p3(input_points[(j + 1) % input_points.size()].x, input_points[(j + 1) % input_points.size()].y);

					CGAL::Object result = CGAL::intersection(Segment_2(p0, p1), Segment_2(p2, p3));

					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
					{
						vecs.push_back(Vector2d(ipoint->x(), ipoint->y()));
					}
				}
			}

			if (distance > 0)
			{
				for (int i = 0; i < vecs.size(); i++)
				{
					double delta_d = Circuit::FindNearestPointPar(vecs[i], input_points);

					if (abs(delta_d - d) > 0.5)
					{
						if (delta_d > d)
						{
							d = d + 1.0;
						}
						else
						{
							delta_d = delta_d + 1.0;
						}
					}

					if (delta_d > d)
					{
						if (delta_d > 1.0)
							delta_d = delta_d - 1.0;

						std::vector<Vector2d>().swap(vecs);
						return delta_d;
					}
				}
			}
			else
			{
				for (int i = 0; i < vecs.size(); i++)
				{
					double delta_d = Circuit::FindNearestPointPar(vecs[i], input_points);

					if (abs(delta_d - d) > 0.5)
					{
						if (delta_d > d)
						{
							d = d + 1.0;
						}
						else
						{
							delta_d = delta_d + 1.0;
						}
					}

					if (delta_d < d)
					{
						if (delta_d > 1.0)
							delta_d = delta_d - 1.0;

						std::vector<Vector2d>().swap(vecs);
						return delta_d;
					}
				}
			}

			std::vector<Vector2d>().swap(vecs);
			return -1.0;
		}

		//points moving on circuit
		static double DeltaDEuclideanDistance(double d, double distance, std::vector<Vector2d> &input_points)
		{
			Vector2d v = Circuit::GetOnePointFromOffset(d, input_points);

			std::vector<Vector2d> vecs;

			int divided_nb = 50;
			for (int i = 0; i < divided_nb; i++)
			{
				Point_2 p0(v[0] + abs(distance)*sin(i * 2 * PI / (double)divided_nb), v[1] + abs(distance)*cos(i * 2 * PI / (double)divided_nb));
				Point_2 p1(v[0] + abs(distance)*sin((i + 1) * 2 * PI / (double)divided_nb), v[1] + abs(distance)*cos((i + 1) * 2 * PI / (double)divided_nb));

				for (int j = 0; j < input_points.size(); j++)
				{
					Point_2 p2(input_points[j].x, input_points[j].y);
					Point_2 p3(input_points[(j + 1) % input_points.size()].x, input_points[(j + 1) % input_points.size()].y);

					CGAL::Object result = CGAL::intersection(Segment_2(p0, p1), Segment_2(p2, p3));

					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
					{
						vecs.push_back(Vector2d(ipoint->x(), ipoint->y()));
					}
				}
			}

			if (distance > 0)
			{
				for (int i = 0; i < vecs.size(); i++)
				{
					double delta_d = Circuit::FindNearestPointPar(vecs[i], input_points);

					if (abs(delta_d - d) > 0.5)
					{
						if (delta_d > d)
						{
							d = d + 1.0;
						}
						else
						{
							delta_d = delta_d + 1.0;
						}
					}

					if (delta_d > d)
					{
						if (delta_d > 1.0)
							delta_d = delta_d - 1.0;

						std::vector<Vector2d>().swap(vecs);
						return delta_d;
					}
				}
			}
			else
			{
				for (int i = 0; i < vecs.size(); i++)
				{
					double delta_d = Circuit::FindNearestPointPar(vecs[i], input_points);

					if (abs(delta_d - d) > 0.5)
					{
						if (delta_d > d)
						{
							d = d + 1.0;
						}
						else
						{
							delta_d = delta_d + 1.0;
						}
					}

					if (delta_d < d)
					{
						if (delta_d > 1.0)
							delta_d = delta_d - 1.0;

						std::vector<Vector2d>().swap(vecs);
						return delta_d;
					}
				}
			}

			std::vector<Vector2d>().swap(vecs);
			return -1.0;
		}

		static double DeltaDistance(double d, Vector2d v0, double distance, std::vector<Vector2d> &input_points)
		{
			Vector2d v = GetOnePointFromOffset(d, input_points);

			int index = Strip::CheckSamePoint(v, input_points);

			if (index >= 0)
			{
				Vector2d v0 = input_points[(index + input_points.size() - 1) % input_points.size()];
				Vector2d v1 = input_points[(index + input_points.size() + 1) % input_points.size()];

				double angle = AnglePI(Vector2d(v0[0] - v[0], v0[1] - v[1]), Vector2d(v1[0] - v[0], v1[1] - v[1]));

				double delta = DeltaDEuclideanDistance(d, distance / sin(angle), input_points);

				if (delta >= 0)
				{
					return delta;
				}
				else
				{
					return DeltaDEuclideanDistance(d, distance, input_points);
				}
			}
			else
			{

				Vector2d v1, v2;
				bool b = false;

				double total_length = GetTotalLength(input_points);
				double length = 0.0;
				for (int i = 0; i < input_points.size(); i++)
				{
					double l = sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));

					if (length / total_length < d&&d < (length + l) / total_length)
					{
						v1 = input_points[i];
						v2 = input_points[(i + 1) % input_points.size()];

						b = true;
						break;
					}
					length += l;
				}
				assert(b);

				double angle = AnglePI(Vector2d(v0[0] - v[0], v0[1] - v[1]), Vector2d(v1[0] - v2[0], v1[1] - v2[1]));

				double delta = DeltaDEuclideanDistance(d, distance / sin(angle), input_points);

				if (delta >= 0)
				{
					return delta;
				}
				else
				{
					return DeltaDEuclideanDistance(d, distance, input_points);
				}
			}
		}

		static double DeltaDistance(double d, double distance, std::vector<Vector2d> &input_points)
		{
			Vector2d v = GetOnePointFromOffset(d, input_points);
			int index = Strip::CheckSamePoint(v, input_points);
			if (index >= 0)
			{
				Vector2d v0 = input_points[(index + input_points.size() - 1) % input_points.size()];
				Vector2d v1 = input_points[(index + input_points.size() + 1) % input_points.size()];

				double angle = AnglePI(Vector2d(v0[0] - v[0], v0[1] - v[1]), Vector2d(v1[0] - v[0], v1[1] - v[1]));
				double delta = DeltaDEuclideanDistance(d, distance / sin(angle), input_points);

				if (delta >= 0)
				{
					return delta;
				}
				else
				{
					return DeltaDEuclideanDistance(d, distance, input_points);
				}

				/*
				double middle_d = d + 0.5;
				if (middle_d > 1.0)
				middle_d = middle_d - 1.0;

				std::vector<Vector2d> part;
				Circuit::SelectOnePartOffset(input_points, d, middle_d, part);

				double start_d = distance;
				for (int i = 0; i < 10; i++)
				{
				start_d = start_d + distance / 10.0;

				double delta = DeltaDEuclideanDistance(d, start_d, input_points);
				Vector2d v_delta = Circuit::GetOnePointFromOffset(delta, input_points);
				double temp_distance = Strip::DistanceLines(v_delta, part);
				if (temp_distance >= distance)
				{
				break;
				}
				}

				return DeltaDEuclideanDistance(d, start_d, input_points);
				*/

			}
			else
			{
				return DeltaDEuclideanDistance(d, distance, input_points);
			}
		}


		static void FindOptimalEntryExitPoints(double toolpath_size, std::vector<Vector2d> &offset, Vector2d &entry_point, Vector2d &exit_point, std::vector<Vector2d> &debug_points)
		{
			double m_d = MAXDOUBLE;
			int m_d_index = -1;

			double m_d_0 = MAXDOUBLE;
			int m_d_index_0 = -1;

			double m_d_1 = MAXDOUBLE;
			int m_d_index_1 = -1;


			for (int i = 0; i < offset.size(); i++)
			{
				Vector2d v0 = offset[(i - 1 + offset.size()) % offset.size()];
				Vector2d v1 = offset[i];
				Vector2d v2 = offset[(i + 1 + offset.size()) % offset.size()];

				double angle = Circuit::Angle2PI(Vector2d(v0[0] - v1[0], v0[1] - v1[1]), Vector2d(v2[0] - v1[0], v2[1] - v1[1]));


				if (abs(angle - PI) < m_d_1)
				{
					m_d_1 = abs(angle - PI);
					m_d_index_1 = i;
				}

				if (angle <= PI)
				{
					if (abs(angle - PI / 2.0) < m_d)
					{
						m_d = abs(angle - PI / 2.0);
						m_d_index = i;
					}
				}

				if (angle <= PI / 2.0)
				{
					if (abs(angle - PI / 2.0) < m_d_0)
					{
						m_d_0 = abs(angle - PI / 2.0);
						m_d_index_0 = i;
					}
				}
			}

			if (m_d_index_0 >= 0 || m_d_index >= 0)
			{
				if (m_d_index_0 >= 0)
				{
					entry_point = offset[m_d_index_0];
				}
				else
				{
					entry_point = offset[m_d_index];
				}
			}
			else
			{
				entry_point = offset[m_d_index_1];
			}

			double d = Circuit::FindNearestPointPar(entry_point, offset);
			//double entry_point_d = Circuit::DeltaDistance(d, toolpath_size, offset);
			
			double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);

			exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);

		}
		static void FindOptimalEntryExitPoints_Richard(double toolpath_size, std::vector<Vector2d> &offset, std::vector<Vector2d> &offset_center, Vector2d &entry_point, Vector2d &exit_point, std::vector<Vector2d> &debug_points)
		{
			Vector2d center(0.0,0.0);

			for (int i = 0; i < offset_center.size(); i++)
			{
				center[0] += offset_center[i][0];
				center[1] += offset_center[i][1];
			}

			center[0] = center[0] / offset_center.size();
			center[1] = center[1] / offset_center.size();

			double par = Circuit::FindNearestPointPar(center, offset);
			entry_point = Circuit::GetOnePointFromOffset(par,offset);

			double d = Circuit::FindNearestPointPar(entry_point, offset);
			//double entry_point_d = Circuit::DeltaDistance(d, toolpath_size, offset);

			double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);

			exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);

		}


		static void FindOptimalEntryExitPoints_Simple(double toolpath_size, std::vector<Vector2d> &offset, std::vector<Vector2d> &next_offset,
			std::vector<Vector2d> &offsets_parts, std::vector<Vector2d> &next_offsets_parts,
			Vector2d &entry_point, Vector2d &exit_point, std::vector<Vector2d> &debug_points)

		{
			for (int i = 0; i < offsets_parts.size(); i++)
			{
				i = offsets_parts.size() / 2.0;

				Vector2d v0 = offsets_parts[i];
				double par = Circuit::FindNearestPointPar(v0, offset);
				double par_exit = DeltaDEuclideanDistance(par, toolpath_size, offset);
				Vector2d v1 = Circuit::GetOnePointFromOffset(par_exit, offset);

				Vector2d connecting_entry_point, connecting_exit_point;
				bool b = Circuit::ComputeNextEntryExitPointForOuter_Richard(toolpath_size, offset, next_offset,
					v0, v1, connecting_entry_point, connecting_exit_point);

				double d1 = Strip::Distance(v1, offsets_parts);


				double d3 = 100.0, d4 = 100.0;
				if (b)
				{
					d3 = Strip::Distance(connecting_entry_point, next_offsets_parts);
					d4 = Strip::Distance(connecting_exit_point, next_offsets_parts);
				}

				if (d1<0.0001&&b&&d3 < 0.0001&&d4 < 0.0001)
				{
					entry_point = v0;
					exit_point = v1;
					return;
				}

				break;
			}
		}
		
		static void FindOptimalEntryExitPoints(double toolpath_size, std::vector<Vector2d> &offset, std::vector<Vector2d> &next_offset,
			std::vector<Vector2d> &offsets_parts, std::vector<Vector2d> &next_offsets_parts,
			Vector2d &entry_point, Vector2d &exit_point, std::vector<Vector2d> &debug_points)
		{
			/*
			double length = Strip::GetTotalLength(offsets_parts);
			
			if (length > toolpath_size)
			{
				double delta = toolpath_size / 2.0 / length;
				
				entry_point = Strip::GetOnePointFromStrip(0.5 + delta, offsets_parts);
				double d = Circuit::FindNearestPointPar(entry_point, offset);
				double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);
				exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);
			}
			else
			{
				entry_point = Strip::GetOnePointFromStrip(0.5, offsets_parts);
				double d = Circuit::FindNearestPointPar(entry_point, offset);
				double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);
				exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);
			}

	
			double d0 = Strip::Distance(entry_point, next_offsets_parts);
			double d1 = Strip::Distance(exit_point, next_offsets_parts);

			if (abs(toolpath_size - d0) > 0.0001 || abs(toolpath_size - d1) > 0.0001)
			{
				length = Strip::GetTotalLength(next_offsets_parts);

				if (length > toolpath_size)
				{
					double delta = toolpath_size / 2.0 / length;
					Vector2d next_entry_point = Strip::GetOnePointFromStrip(0.5, next_offsets_parts);

					debug_points.push_back(next_entry_point);

					double d = Circuit::FindNearestPointPar(next_entry_point, offset);
					double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);
					entry_point = Circuit::GetOnePointFromOffset(d, offset);
					exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);
				}
				else
				{
					Vector2d next_entry_point = Strip::GetOnePointFromStrip(0.5, next_offsets_parts);
					double d = Circuit::FindNearestPointPar(next_entry_point, offset);
					double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);
					entry_point = Circuit::GetOnePointFromOffset(d, offset);
					exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);
				}

			}
			*/

			double m_d = MAXDOUBLE;
			int m_d_index = -1;

			double m_d_0 = MAXDOUBLE;
			int m_d_index_0 = -1;

			double m_d_1 = MAXDOUBLE;
			int m_d_index_1 = -1;

			for (int i = 0; i < offset.size(); i++)
			{
				double d1 = Strip::Distance(offset[i], offsets_parts);

				double d = Circuit::FindNearestPointPar(offset[i], offset);
				double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);
				Vector2d v_exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);

				double d2 = Strip::Distance(v_exit_point, offsets_parts);

				Vector2d connecting_entry_point, connecting_exit_point;

				bool b=Circuit::ComputeNextEntryExitPointForOuter_Richard(toolpath_size, offset, next_offset,
					offset[i], v_exit_point, connecting_entry_point, connecting_exit_point);

				double d3=100.0, d4=100.0;
				if (b)
				{
					d3 = Strip::Distance(connecting_entry_point, next_offsets_parts);
					d4 = Strip::Distance(connecting_exit_point, next_offsets_parts);
				}

				if (d1<0.0001&&d2<0.0001&&b&&d3<0.0001&&d4<0.0001)
				{
					Vector2d v0 = offset[(i - 1 + offset.size()) % offset.size()];
					Vector2d v1 = offset[i];
					Vector2d v2 = offset[(i + 1 + offset.size()) % offset.size()];

					double angle = Circuit::Angle2PI(Vector2d(v0[0] - v1[0], v0[1] - v1[1]), Vector2d(v2[0] - v1[0], v2[1] - v1[1]));

					if (abs(angle - PI) < m_d_1)
					{
						m_d_1 = abs(angle - PI);
						m_d_index_1 = i;
					}

					if (angle <= PI)
					{
						if (abs(angle - PI / 2.0) < m_d)
						{
							m_d = abs(angle - PI / 2.0);
							m_d_index = i;
						}
					}

					if (angle <= PI / 2.0)
					{
						if (abs(angle - PI / 2.0) < m_d_0)
						{
							m_d_0 = abs(angle - PI / 2.0);
							m_d_index_0 = i;
						}
					}
				}
			}

			if (m_d_index_0 >= 0 || m_d_index >= 0)
			{
				if (m_d_index_0 >= 0)
				{
					entry_point = offset[m_d_index_0];
				}
				else
				{
					entry_point = offset[m_d_index];
				}
			}
			else
			{
				entry_point = offset[m_d_index_1];
			}

			double d = Circuit::FindNearestPointPar(entry_point, offset);
			double entry_point_d = DeltaDEuclideanDistance(d, toolpath_size, offset);
			exit_point = Circuit::GetOnePointFromOffset(entry_point_d, offset);
		}


	static Vector2d DetectNearestPoint(Vector2d v0, Vector2d v1, std::vector<Vector2d> &input_points)
			{
			std::vector<Vector2d> vecs;

			for (int i = 0; i < input_points.size(); i++)
			{
				Vector2d v2 = input_points[i];
				Vector2d v3 = input_points[(i + 1) % input_points.size()];

				CGAL::Object result = CGAL::intersection(Line_2(Point_2(v0[0], v0[1]), Point_2(v1[0], v1[1])), Segment_2(Point_2(v2[0], v2[1]), Point_2(v3[0], v3[1])));

				if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
				{
					vecs.push_back(Vector2d(ipoint->x(), ipoint->y()));
				}
			}

			assert(vecs.size() != 0);

			double min_d = Strip::Distance(v0, vecs[0]);
			int index = 0;

			for (int i = 1; i < vecs.size(); i++)
			{
				double d = Strip::Distance(v0, vecs[i]);

				if (d < min_d)
				{
					min_d = d;
					index = i;
				}
			}

			return vecs[index];
		}


	static bool DetectNearestPoint(Vector2d v0, Vector2d v1, std::vector<Vector2d> &input_points, Vector2d &output_point)
	{
		std::vector<Vector2d> vecs;

		for (int i = 0; i < input_points.size(); i++)
		{
			Vector2d v2 = input_points[i];
			Vector2d v3 = input_points[(i + 1) % input_points.size()];

			CGAL::Object result = CGAL::intersection(Ray_2(Point_2(v0[0], v0[1]), Point_2(v1[0], v1[1])), Segment_2(Point_2(v2[0], v2[1]), Point_2(v3[0], v3[1])));

			if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
			{
				vecs.push_back(Vector2d(ipoint->x(), ipoint->y()));
			}
		}

		if (vecs.size() == 0)
			return false;


		double min_d = Strip::Distance(v0, vecs[0]);
		int index = 0;

		for (int i = 1; i < vecs.size(); i++)
		{
			double d = Strip::Distance(v0, vecs[i]);

			if (d < min_d)
			{
				min_d = d;
				index = i;
			}
		}

		output_point=vecs[index];

		return true;
	}



		static void ComputeNextEntryExitPointForOuter(double toolpath_size, std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset,
			Vector2d &entry_point, Vector2d &exit_point, Vector2d &next_entry_point, Vector2d &next_exit_point)
		{
			/*
			Vector2d n(exit_point[0] - entry_point[0], exit_point[1] - entry_point[1]);
			Vector2d pn(-n[1], n[0]);

			next_entry_point = DetectNearestPoint(entry_point, Vector2d(entry_point[0] + pn[0], entry_point[1] + pn[1]), inside_offset);
			next_exit_point = DetectNearestPoint(exit_point, Vector2d(exit_point[0] + pn[0], exit_point[1] + pn[1]), inside_offset);
			*/

			next_entry_point=Circuit::FindNearestPoint(entry_point, inside_offset);
			next_exit_point = Circuit::FindNearestPoint(exit_point, inside_offset);

		}


		static bool ComputeNextEntryExitPointForOuter_Richard(double toolpath_size, std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset,
			Vector2d &entry_point, Vector2d &exit_point, Vector2d &next_entry_point, Vector2d &next_exit_point)
		{
			Vector2d n(entry_point[0] - exit_point[0], entry_point[1] - exit_point[1]);
			Vector2d pn(-n[1], n[0]);

			bool b0 = DetectNearestPoint(entry_point, Vector2d(entry_point[0] + pn[0], entry_point[1] + pn[1]), inside_offset, next_entry_point);
			bool b1 = DetectNearestPoint(exit_point, Vector2d(exit_point[0] + pn[0], exit_point[1] + pn[1]), inside_offset, next_exit_point);

			if (b0&&b1)
			{
				return true;
			}
			else
			{
				return false;
			}
		}


		//	static void ComputeNextEntryExitPointForOuter(std::vector<Vector2d> &input_points, Vector2d &)

		static Vector2d FindNearestPoint(Vector2d &v, std::vector<Vector2d> &input_points)
		{
			double d = Circuit::FindNearestPointPar(v,input_points);
			return Circuit::GetOnePointFromOffset(d,input_points);
		}

		static void FindNearestPoints(Vector2d &v, std::vector<Vector2d> &input_points, std::vector<Vector2d> &nearest_points)
		{
			double distance = Circuit::Distance(v, input_points);
		
			for (int i = 0; i < input_points.size(); i++)
			{
				Vector2d v0 = input_points[i];
				Vector2d v1 = input_points[(i + 1) % input_points.size()];

				double d = Strip::Distance(v, v0, v1);

				if (abs(d - distance) < 0.0001)
				{
					nearest_points.push_back(Strip::FindNearestPoint(v, v0, v1));
				}
			}
		}

		static Vector2d FindNearestPoint(std::vector<Vector2d> &input_points, Vector2d &ray_start, Vector2d &ray_end)
		{
			std::vector<Vector2d> vecs;
			for (int i = 0; i < input_points.size(); i++)
			{
				Vector2d v0 = input_points[i];
				Vector2d v1 = input_points[(i + 1) % input_points.size()];

				CGAL::Object result = CGAL::intersection(Segment_2(Point_2(v0[0], v0[1]), Point_2(v1[0], v1[1])), Line_2(Point_2(ray_start[0], ray_start[1]), Point_2(ray_end[0], ray_end[1])));

				if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
				{
					vecs.push_back(Vector2d(ipoint->x(), ipoint->y()));
				}
			}

			assert(vecs.size() != 0);

			double min_d = Strip::Distance(vecs[0], ray_end);
			int min_index = 0;

			for (int i = 1; i < vecs.size(); i++)
			{
				double d = Strip::Distance(vecs[i], ray_end);
				if (d < min_d)
				{
					min_d = d;
					min_index = i;
				}
			}
			Vector2d v = vecs[min_index];
			std::vector<Vector2d>().swap(vecs);
			return v;
		}

		static void ComputeNextEntryExitPointForInner(double toolpath_size, std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset,
			Vector2d &entry_point, Vector2d &exit_point, Vector2d &next_entry_point, Vector2d &next_exit_point)
		{
			double entry_d_1 = Circuit::FindNearestPointPar(entry_point, inside_offset);
			double exit_d_1 = Circuit::FindNearestPointPar(exit_point, inside_offset);

			Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, inside_offset);
			Vector2d exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, inside_offset);

			//handle the case where the next_exit_d_0 and next_entry_d_0 meets each other

			int index_0 = Strip::CheckSamePoint(entry_p_1, inside_offset);
			int index_1 = Strip::CheckSamePoint(exit_p_1, inside_offset);

			//if (abs(entry_d_1 - exit_d_1) < 0.00001)
			if (index_0 >= 0 || index_1 >= 0)
			{
				Vector2d n_0(entry_point[0] - entry_p_1[0], entry_point[1] - entry_p_1[1]);
				Vector2d n_1(exit_point[0] - entry_p_1[0], exit_point[1] - entry_p_1[1]);

				if (n_0[0] * n_1[1] - n_0[1] * n_1[0]>0)
				{
					if (index_0 >= 0)
					{
						double temp = Circuit::DeltaDistance(entry_d_1, toolpath_size, inside_offset);
						if (temp >= 0)
							exit_d_1 = temp;
					}
					else
					{
						//double temp = Circuit::DeltaDistance(exit_d_1, -toolpath_size, inside_offset);
						//if (temp >= 0)
						//	entry_d_1 = temp;
					}
				}
				else
				{
					if (index_1 >= 0)
					{
						double temp = Circuit::DeltaDistance(exit_d_1, toolpath_size, inside_offset);
						if (temp >= 0)
							entry_d_1 = temp;
					}
					else
					{
						//	double temp = Circuit::DeltaDistance(entry_d_1, -toolpath_size, inside_offset);
						//	if (temp >= 0)
						//	exit_d_1 = temp;
					}

				}
				entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, inside_offset);
				exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, inside_offset);
			}

			next_entry_point = entry_p_1;
			next_exit_point = exit_p_1;
		}

		static void ComputeNextEntryExitPointForInner_Richard(double toolpath_size, std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset,
			Vector2d &entry_point, Vector2d &exit_point, Vector2d &next_entry_point, Vector2d &next_exit_point)
		{
			double entry_d_1 = Circuit::FindNearestPointPar(entry_point, inside_offset);
			double exit_d_1 = Circuit::FindNearestPointPar(exit_point, inside_offset);

			Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, inside_offset);
			Vector2d exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, inside_offset);

			//handle the case where the next_exit_d_0 and next_entry_d_0 meets each other

			int index_0 = Strip::CheckSamePoint(entry_p_1, inside_offset);
			int index_1 = Strip::CheckSamePoint(exit_p_1, inside_offset);

			//if (abs(entry_d_1 - exit_d_1) < 0.00001)
			if (index_0 >= 0 || index_1 >= 0)
			{
				Vector2d n_0(entry_point[0] - entry_p_1[0], entry_point[1] - entry_p_1[1]);
				Vector2d n_1(exit_point[0] - entry_p_1[0], exit_point[1] - entry_p_1[1]);

				if (n_0[0] * n_1[1] - n_0[1] * n_1[0]>0)
				{
					if (index_0 >= 0)
					{
						//double temp = Circuit::DeltaDistance(entry_d_1, toolpath_size, inside_offset);
						double temp = DeltaDEuclideanDistance(entry_d_1, toolpath_size, inside_offset);
						if (temp >= 0)
							exit_d_1 = temp;
					}
					else
					{
						//double temp = Circuit::DeltaDistance(exit_d_1, -toolpath_size, inside_offset);
						//if (temp >= 0)
						//	entry_d_1 = temp;
					}
				}
				else
				{
					if (index_1 >= 0)
					{
						//double temp = Circuit::DeltaDistance(exit_d_1, toolpath_size, inside_offset);
						double temp = DeltaDEuclideanDistance(exit_d_1, toolpath_size, inside_offset);
						if (temp >= 0)
							entry_d_1 = temp;
					}
					else
					{
						//	double temp = Circuit::DeltaDistance(entry_d_1, -toolpath_size, inside_offset);
						//	if (temp >= 0)
						//	exit_d_1 = temp;
					}

				}
				entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, inside_offset);
				exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, inside_offset);
			}

			next_entry_point = entry_p_1;
			next_exit_point = exit_p_1;
		}



		static void ComputeNextEntryExitPointForInner_Update_Richard(int index, double toolpath_size, std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset,
			Vector2d &entry_point, Vector2d &exit_point, Vector2d &update_entry_point, Vector2d &update_exit_point)
		{
			double entry_d_1 = Circuit::FindNearestPointPar(entry_point, inside_offset);
			double exit_d_1 = Circuit::FindNearestPointPar(exit_point, inside_offset);

			Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, inside_offset);
			Vector2d exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, inside_offset);

			update_entry_point = entry_point;
			update_exit_point = exit_point;

			//handle the case where the next_exit_d_0 and next_entry_d_0 meets each other
			int index_0 = Strip::CheckSamePoint(entry_p_1, inside_offset);
			int index_1 = Strip::CheckSamePoint(exit_p_1, inside_offset);

			//if (abs(entry_d_1 - exit_d_1) < 0.00001)
			if (index_0 >= 0 && index_1 >= 0 && index_0 == index_1)
			{
				double temp = DeltaDEuclideanDistance(exit_d_1, toolpath_size, inside_offset);

				if (temp >= 0)
				{
					double entry_d_0 = Circuit::FindNearestPointPar(entry_point, outside_offset);
					double exit_d_0 = Circuit::FindNearestPointPar(exit_point, outside_offset);

					std::vector<Vector2d> one_path_0;
					std::vector<Vector2d> one_path_1;
					Circuit::SelectOnePartOffset(outside_offset, entry_d_0, exit_d_0, one_path_0);
					Circuit::SelectOnePartOffset(outside_offset, exit_d_0, entry_d_0, one_path_1);

					if (Strip::GetTotalLength(one_path_0) < Strip::GetTotalLength(one_path_1))
					{
						//entry
						Vector2d p = Circuit::GetOnePointFromOffset(temp, inside_offset);
						double par = Circuit::FindNearestPointPar(p, outside_offset);
						update_entry_point = Circuit::GetOnePointFromOffset(par, outside_offset);
					}
					else
					{
						//exit
						Vector2d p = Circuit::GetOnePointFromOffset(temp, inside_offset);
						double par = Circuit::FindNearestPointPar(p, outside_offset);
						update_exit_point = Circuit::GetOnePointFromOffset(par, outside_offset);
					}

				}

			}
		}



		//for fermat spiral generation
		static Line_2 GetRelatedLine(Vector2d v, std::vector<Vector2d> &input_points)
		{
			assert(Distance(v, input_points)<0.0001);

			Line_2 line_2;

			int index = Strip::CheckSamePoint(v, input_points);

			if (index >= 0)
			{
				Vector2d v1 = input_points[(index + input_points.size() - 1) % input_points.size()];
				line_2 = Line_2(Point_2(v[0], v[1]), Point_2(v1[0], v1[1]));
			}
			else
			{
				double d = GetTotalLength(input_points);
				double par = FindNearestPointPar(v, input_points);
				double length = 0.0;
				for (int i = 0; i < input_points.size(); i++)
				{
					double l = sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));

					if (length / d < par&&par < (length + l) / d)
					{
						Vector2d v1 = input_points[i];
						Vector2d v2 = input_points[(i + input_points.size() - 1) % input_points.size()];

						Vector2d v3 = v;
						v3[0] = v3[0] + v2[0] - v1[0];
						v3[1] = v3[1] + v2[1] - v1[1];

						line_2 = Line_2(Point_2(v[0], v[1]), Point_2(v3[0], v3[1]));

						break;
					}
					length += l;
				}
			}
			return line_2;
		}

		static Line_2 GetTangent(Vector2d v, std::vector<Vector2d> &input_points)
		{
			assert(Distance(v, input_points)<0.0001);

			int index = Strip::CheckSamePoint(v, input_points);

			Line_2 line_2;

			if (index >= 0)
			{
				Vector2d v0 = input_points[(index + input_points.size() - 1) % input_points.size()];
				Vector2d v1 = input_points[(index + input_points.size() + 1) % input_points.size()];

				Vector2d n(v0[0] - v[0], v0[1] - v[1]);
				double n_l = std::sqrt(n[0] * n[0] + n[1] * n[1]);
				v0[0] = v[0] + n[0] / n_l;
				v0[1] = v[1] + n[1] / n_l;

				n = Vector2d(v1[0] - v[0], v1[1] - v[1]);
				n_l = std::sqrt(n[0] * n[0] + n[1] * n[1]);
				v1[0] = v[0] + n[0] / n_l;
				v1[1] = v[1] + n[1] / n_l;

				line_2 = Line_2(Point_2(v[0] + v1[0] - v0[0], v[1] + v1[1] - v0[1]), Point_2(v[0], v[1]));
			}
			else
			{
				double d = GetTotalLength(input_points);
				double par = FindNearestPointPar(v, input_points);
				double length = 0.0;
				for (int i = 0; i < input_points.size(); i++)
				{
					double l = sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
					if (length / d < par&&par < (length + l) / d)
					{
						Vector2d v0 = input_points[i];
						Vector2d v1 = input_points[(i + 1) % input_points.size()];

						line_2 = Line_2(Point_2(v0[0], v0[1]), Point_2(v1[0], v1[1]));
						break;
					}
					length += l;
				}
			}

			return line_2;
		}

		//for offset graph decomposition
		static bool compare_double(const double &d0, const double &d1)
		{
			return d0 > d1;
		}


		static void CutACircuit(std::vector<Vector2d> &input_points, TrunkNode &trunk_node)
		{
			std::vector<double> cutting_points_par;

			for (int i = 0; i <trunk_node.connecting_points.size(); i = i + 2)
			{
				//double par = (trunk_node.connecting_points[i] + trunk_node.connecting_points[i + 1]) / 2.0;

				Vector2d v0 = Circuit::GetOnePointFromOffset(trunk_node.connecting_points[i], input_points);
				Vector2d v1 = Circuit::GetOnePointFromOffset(trunk_node.connecting_points[i + 1], input_points);
				Vector2d v3((v0[0] + v1[0]) / 2.0, (v0[1] + v1[1]) / 2.0);
				double par = Circuit::FindNearestPointPar(v3, input_points);

				cutting_points_par.push_back(par);
			}

			std::sort(cutting_points_par.begin(), cutting_points_par.end(), compare_double);

			if (cutting_points_par.size() > 1)
			{
				for (int i = 0; i < cutting_points_par.size(); i++)
				{
					double cur_par = cutting_points_par[i];
					double next_par = cutting_points_par[(i + 1) % cutting_points_par.size()];

					int cur_par_index = -1;
					int next_par_index = -1;

					for (int j = 0; j < trunk_node.connecting_points.size() && (cur_par_index<0 || next_par_index<0); j = j + 2)
					{
						//double par = (trunk_node.connecting_points[j] + trunk_node.connecting_points[j + 1]) / 2.0;

						Vector2d v0 = Circuit::GetOnePointFromOffset(trunk_node.connecting_points[j], input_points);
						Vector2d v1 = Circuit::GetOnePointFromOffset(trunk_node.connecting_points[j + 1], input_points);
						Vector2d v3((v0[0] + v1[0]) / 2.0, (v0[1] + v1[1]) / 2.0);
						double par = Circuit::FindNearestPointPar(v3, input_points);


						if (abs(par - cur_par) < 0.00001)
						{
							cur_par_index = j;
						}

						if (abs(par - next_par) < 0.00001)
						{
							next_par_index = j;
						}
					}

					if (cur_par_index >= 0 && next_par_index >= 0)
					{
						trunk_node.connecting_points_pairs.push_back(cur_par_index + 1);
						trunk_node.connecting_points_pairs.push_back(next_par_index);
					}

				}

			}
			else
			{
				if (cutting_points_par.size() == 1)
				{
					trunk_node.connecting_points_pairs.push_back(0);
					trunk_node.connecting_points_pairs.push_back(1);
				}
			}
		
			std::vector<double>().swap(cutting_points_par);
		}


		static void CutACircuit(std::vector<Vector2d> &input_points, std::vector<double> &cutting_points,
			std::vector<std::vector<Vector2d>> &pathes, std::vector<int> &connecting)
		{
			std::vector<double> cutting_points_par;

			for (int i = 0; i < cutting_points.size(); i = i + 2)
			{
				double par = (cutting_points[i] + cutting_points[i + 1]) / 2.0;
				cutting_points_par.push_back(par);
			}

			std::sort(cutting_points_par.begin(), cutting_points_par.end(), compare_double);


			for (int i = 0; i < cutting_points_par.size(); i++)
			{
				double cur_par = cutting_points_par[i];
				double next_par = cutting_points_par[(i + 1) % cutting_points_par.size()];

				int cur_par_index = -1;
				int next_par_index = -1;

				for (int j = 0; j < cutting_points.size() && (cur_par_index<0 || next_par_index<0); j = j + 2)
				{
					double par = (cutting_points[j] + cutting_points[j + 1]) / 2.0;

					if (abs(par - cur_par) < 0.00001)
					{
						cur_par_index = j;
					}

					if (abs(par - next_par) < 0.00001)
					{
						next_par_index = j;
					}
				}

				if (cur_par_index >= 0 && next_par_index >= 0)
				{
					//SelectOnePartOffset();
					std::vector<Vector2d> one_path;

					//SelectOnePartOffset(input_points, cutting_points[cur_par_index + 1], cutting_points[next_par_index], one_path);

					connecting.push_back(cur_par_index + 1);
					connecting.push_back(next_par_index);

					pathes.push_back(one_path);

					std::vector<Vector2d>().swap(one_path);
				}

			}
			std::vector<double>().swap(cutting_points_par);
		}


		static void DetectSharingPart0(double toolpath_size, std::vector<Vector2d> &offset_0, std::vector<Vector2d> &offset_1,
			std::vector<std::vector<Vector2d>> &offset_0_path, std::vector<std::vector<Vector2d>> &offset_1_path)
		{
			bool b0 = CheckInside(offset_0[0], offset_1);
			bool b1 = CheckInside(offset_1[0], offset_0);

			if (b0)
			{
				GenerateOffsetForOutsideWithClipper(offset_0, toolpath_size, offset_1_path);
				GenerateOffsetWithClipper(offset_1, toolpath_size, offset_0_path);
			}

			if (b1)
			{
				GenerateOffsetForOutsideWithClipper(offset_1, toolpath_size, offset_0_path);
				GenerateOffsetWithClipper(offset_0, toolpath_size, offset_1_path);
			}

			if (!b0&&!b1)
			{
				GenerateOffsetForOutsideWithClipper(offset_0, toolpath_size, offset_1_path);
				GenerateOffsetForOutsideWithClipper(offset_1, toolpath_size, offset_0_path);
			}
		}

		static void  OffsetRelatedContour(double offset_distance, std::vector<Vector2d> &offset, Polygon_with_holes &contours,
			std::vector<std::vector<Vector2d>> &offset_parts,
			std::vector<Vector2d> &debug_points, std::vector<std::vector<Vector2d>> &debug_lines)
		{

			//get contour points
			std::vector<std::vector<Vector2d>> contour_points;
			std::vector<Vector2d> contour;
			for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
			{
				contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
			}
			contour_points.push_back(contour);

			for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
			{
				std::vector<Vector2d>().swap(contour);
				for (Polygon_2::Vertex_iterator ver_iter = hole_iter->vertices_begin(); ver_iter != hole_iter->vertices_end(); ver_iter++)
				{
					contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				}
				contour_points.push_back(contour);
			}
			std::vector<Vector2d>().swap(contour);


			//initial offset_related_contour_index
			std::vector<std::vector<int>> offset_related_contour_index(offset.size(), std::vector<int>());
			std::vector<std::vector<Vector2d>> offset_related_contour_point(offset.size(), std::vector<Vector2d>());
			std::vector<bool>  offset_related_contour_used(offset.size(), true);
			std::vector<int> cutting_point_index;

			for (int i = 0; i < offset.size(); i++)
			{
				for (int j = 0; j < contour_points.size(); j++)
				{
					double d0 = Circuit::Distance(offset[i], contour_points[j]);
					double d1 = abs(d0 - offset_distance);
					if (d1<0.01)
					{
						std::vector<Vector2d> ps;

						Circuit::FindNearestPoints(offset[i], contour_points[j], ps);

						for (int k = 0; k < ps.size(); k++)
						{
							offset_related_contour_point[i].push_back(ps[k]);
							offset_related_contour_index[i].push_back(j);
						}

						std::vector<Vector2d>().swap(ps);
					}
				}

				if (offset_related_contour_index[i].size() >= 2)
				{
					cutting_point_index.push_back(i);

					for (int j = 0; j < offset_related_contour_point[i].size(); j++)
					{
					}
				}
			}

			//get the related parts
			std::vector<std::vector<double>> offset_related_contour_parts(contour_points.size(), std::vector<double>());

			if (cutting_point_index.size()==0)
			{
				offset_related_contour_parts[0].push_back(0.0);
				offset_related_contour_parts[0].push_back(1.0);
			}
			else
			{
				for (int i = 0; i < cutting_point_index.size(); i++)
				{
					int start_index = cutting_point_index[i];
					int next_index = cutting_point_index[(i + 1) % cutting_point_index.size()];

					if (next_index != (start_index + 1) % offset.size())
					{
						int related_index = offset_related_contour_index[(start_index + 1) % offset.size()][0];

						Vector2d start_index_point, next_index_point;

						for (int j = 0; j < offset_related_contour_index[start_index].size(); j++)
						{
							if (offset_related_contour_index[start_index][j] == related_index)
							{
								start_index_point = offset_related_contour_point[start_index][j];
							}
						}

						for (int j = 0; j < offset_related_contour_index[next_index].size(); j++)
						{
							if (offset_related_contour_index[next_index][j] == related_index)
							{
								next_index_point = offset_related_contour_point[next_index][j];
							}
						}
					
						offset_related_contour_parts[related_index].push_back(Circuit::FindNearestPointPar(next_index_point, contour_points[related_index]));
						offset_related_contour_parts[related_index].push_back(Circuit::FindNearestPointPar(start_index_point, contour_points[related_index]));
					}
				}
			}

			for (int i = 0; i < offset_related_contour_parts.size(); i++)
			{
				for (int j = 0; j < offset_related_contour_parts[i].size(); j = j + 2)
				{
					if (abs(offset_related_contour_parts[i][j]) < 0.00001&&abs(offset_related_contour_parts[i][j + 1]-1.0) < 0.00001)
					{
						offset_parts.push_back(contour_points[i]);
					}
					else
					{
						std::vector<Vector2d> one_path;
						Circuit::SelectOnePartOffset(contour_points[i], offset_related_contour_parts[i][j], offset_related_contour_parts[i][j + 1], one_path);
						offset_parts.push_back(one_path);
					}
				
				}
			}

			std::vector<std::vector<double>>().swap(offset_related_contour_parts);
			std::vector<int>().swap(cutting_point_index);
			std::vector<std::vector<int>>().swap(offset_related_contour_index);
			std::vector<std::vector<Vector2d>>().swap(offset_related_contour_point);
			std::vector<bool>().swap(offset_related_contour_used);
			std::vector<std::vector<Vector2d>>().swap(contour_points);
		}

		static bool PointsOnSharingEdge(Vector2d v, std::vector<std::vector<Vector2d>> &offsets_parts_0,
			std::vector<std::vector<Vector2d>> &offsets_parts_1 )
		{
			Vector2d p0 = Strip::FindNearestPoint(v, offsets_parts_0);
			Vector2d p1 = Strip::FindNearestPoint(p0, offsets_parts_1);
			double d = Strip::Distance(p1,p0);

			return abs(d) < 0.001;
		}

		static double Small_cut_direction(std::vector<Vector2d> &input_points, double d0, double d1)
		{
			std::vector<Vector2d> one_path;
			std::vector<Vector2d> another_path;

			Circuit::SelectOnePartOffset(input_points, d0, d1, one_path);
			Circuit::SelectOnePartOffset(input_points, d1, d0, another_path);

			double length_0 = Strip::GetTotalLength(one_path);
			double length_1 = Strip::GetTotalLength(another_path);

			if (length_0 < length_1)
			{
				return d1;
			}
			else
			{
				return d0;
			}

		}

		static bool ConnectTwoTrunkNodes_Richard(double toolpath_size, std::vector<Vector2d> &offset0, std::vector<Vector2d> &offset1,
			int trunk_node_index_0, int trunk_node_index_1,
			TrunkNode &trunk_node_0, TrunkNode &trunk_node_1,
			std::vector<Vector2d> &offsets_parts_0,
			std::vector<Vector2d> &offsets_parts_1,
			std::vector<Vector2d> &debug_points, std::vector<std::vector<Vector2d>> &debug_lines, int index)
		{
			bool b = false;
			for (int i = 0; i < offsets_parts_0.size() && !b; i++)
			{
				Vector2d current_p = offsets_parts_0[i];

				if (abs(Circuit::Distance(current_p, offset1) - toolpath_size)>0.0001)
					continue;

				double current_par = Circuit::FindNearestPointPar(current_p, offset0);
				double neighbor_par = Circuit::DeltaDEuclideanDistance(current_par, toolpath_size, offset0);
				Vector2d neighbor_p = Circuit::GetOnePointFromOffset(neighbor_par, offset0);

				if (abs(Circuit::Distance(neighbor_p, offset1) - toolpath_size)>0.0001)
					continue;

				if (Strip::Distance(neighbor_p, offsets_parts_0) > 0.00001)
				{
					continue;
				}

				double rel_current_par = Circuit::FindNearestPointPar(current_p, offset1);
				double rel_neighbor_par = Circuit::FindNearestPointPar(neighbor_p, offset1);

				Vector2d rel_current_p = Circuit::GetOnePointFromOffset(rel_current_par, offset1);
				Vector2d rel_neighbor_p = Circuit::GetOnePointFromOffset(rel_neighbor_par, offset1);

				if (abs(Circuit::Distance(rel_current_p, offset0) - toolpath_size)>0.0001)
					continue;

				if (abs(Circuit::Distance(rel_neighbor_p, offset0) - toolpath_size)>0.0001)
					continue;

				if (Strip::Distance(rel_current_p, offsets_parts_1) < 0.00001&&
					Strip::Distance(rel_neighbor_p, offsets_parts_1) < 0.00001)
				{
					trunk_node_0.connecting_trunk_nodes_points.push_back(current_par);
					trunk_node_0.connecting_trunk_nodes_points.push_back(neighbor_par);

					trunk_node_1.connecting_trunk_nodes_points.push_back(rel_current_par);
					trunk_node_1.connecting_trunk_nodes_points.push_back(rel_neighbor_par);

					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 2);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 1);

					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 2);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 1);

					b = true;
				}
			}

			for (int i = 0; i < offsets_parts_1.size() && !b; i++)
			{
				Vector2d current_p = offsets_parts_1[i];

				if (abs(Circuit::Distance(current_p, offset0) - toolpath_size)>0.0001)
					continue;

				double current_par = Circuit::FindNearestPointPar(current_p, offset1);
				double neighbor_par = Circuit::DeltaDEuclideanDistance(current_par, toolpath_size, offset1);
				Vector2d neighbor_p = Circuit::GetOnePointFromOffset(neighbor_par, offset1);


				if (abs(Circuit::Distance(neighbor_p, offset0) - toolpath_size)>0.0001)
					continue;

				if (Strip::Distance(neighbor_p, offsets_parts_1) > 0.00001)
				{
					continue;
				}

				double rel_current_par = Circuit::FindNearestPointPar(current_p, offset0);
				double rel_neighbor_par = Circuit::FindNearestPointPar(neighbor_p, offset0);

				Vector2d rel_current_p = Circuit::GetOnePointFromOffset(rel_current_par, offset0);
				Vector2d rel_neighbor_p = Circuit::GetOnePointFromOffset(rel_neighbor_par, offset0);

				if (abs(Circuit::Distance(rel_current_p, offset1) - toolpath_size)>0.0001)
					continue;

				if (abs(Circuit::Distance(rel_neighbor_p, offset1) - toolpath_size)>0.0001)
					continue;

				if (Strip::Distance(rel_current_p, offsets_parts_0) < 0.00001&&
					Strip::Distance(rel_neighbor_p, offsets_parts_0) < 0.00001)
				{
					trunk_node_1.connecting_trunk_nodes_points.push_back(neighbor_par);
					trunk_node_1.connecting_trunk_nodes_points.push_back(current_par);
		
					trunk_node_0.connecting_trunk_nodes_points.push_back(rel_neighbor_par);
					trunk_node_0.connecting_trunk_nodes_points.push_back(rel_current_par);

					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 2);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 1);

					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 2);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 1);

					b = true;
				}
			}

			return b;
		}


		static bool ConnectTwoTrunkNodes(double toolpath_size, std::vector<Vector2d> &offset0, std::vector<Vector2d> &offset1,
			int trunk_node_index_0, int trunk_node_index_1,
			TrunkNode &trunk_node_0, TrunkNode &trunk_node_1,
			std::vector<Vector2d> &offsets_parts_0,
			std::vector<Vector2d> &offsets_parts_1,
			std::vector<Vector2d> &debug_points, std::vector<std::vector<Vector2d>> &debug_lines,int index)
		{
			//std::vector<std::vector<Vector2d>> offset_0_pathes;
			//std::vector<std::vector<Vector2d>> offset_1_pathes;
		//	DetectSharingPart0(toolpath_size, offset0, offset1, offset_0_pathes, offset_1_pathes);


			bool b = false;
			for (int i = 0; i < trunk_node_0.connecting_leaf_nodes_points.size() && !b; i = i + 2)
			{
				double cp_0 = trunk_node_0.connecting_leaf_nodes_points[i];
				double cp_1 = trunk_node_0.connecting_leaf_nodes_points[i + 1];

				Vector2d v0 = GetOnePointFromOffset(cp_0, offset0);
				Vector2d v1 = GetOnePointFromOffset(cp_1, offset0);

				double d0 = Strip::Distance(v0, offsets_parts_0);
				double d1 = Strip::Distance(v1, offsets_parts_0);

				//bool b0 = PointsOnSharingEdge(v0, offsets_parts_0, offsets_parts_1);
				//bool b1 = PointsOnSharingEdge(v1, offsets_parts_0, offsets_parts_1);

				if (d0 < 0.00001&&d1 < 0.00001)
				//if (b0&&b1)
				{
					double d00 = DeltaDEuclideanDistance(cp_0, toolpath_size, offset0);
					Vector2d v00 = GetOnePointFromOffset(d00, offset0);

					Vector2d next_entry_point;
					Vector2d next_exit_point;

					ComputeNextEntryExitPointForInner(toolpath_size, offset0, offset1, v0, v00, next_entry_point, next_exit_point);

					b = true;

					trunk_node_0.connecting_trunk_nodes_points.push_back(d00);
					trunk_node_0.connecting_trunk_nodes_points.push_back(cp_0);

					trunk_node_1.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_exit_point, offset1));
					trunk_node_1.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_entry_point, offset1));

					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 2);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 1);

					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 2);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 1);

				}
			}

			for (int i = 0; i < trunk_node_1.connecting_leaf_nodes_points.size() && !b; i = i + 2)
			{
				double cp_0 = trunk_node_1.connecting_leaf_nodes_points[i];
				double cp_1 = trunk_node_1.connecting_leaf_nodes_points[i + 1];
				Vector2d v0 = GetOnePointFromOffset(cp_0, offset1);
				Vector2d v1 = GetOnePointFromOffset(cp_1, offset1);

				double d0 = Strip::Distance(v0, offsets_parts_1);
				double d1 = Strip::Distance(v1, offsets_parts_1);

				//bool b0 = PointsOnSharingEdge(v0, offsets_parts_1, offsets_parts_0);
				//bool b1 = PointsOnSharingEdge(v1, offsets_parts_1, offsets_parts_0);

		

				//if (b0&&b1)
				if (d0 < 0.00001&&d1 < 0.00001)
				{

					double d11 = DeltaDEuclideanDistance(cp_1, -toolpath_size, offset1);
					Vector2d v11 = GetOnePointFromOffset(d11, offset1);


					Vector2d next_entry_point;
					Vector2d next_exit_point;

					ComputeNextEntryExitPointForOuter(toolpath_size, offset1, offset0, v1, v11, next_entry_point, next_exit_point);

		

					b = true;

					trunk_node_1.connecting_trunk_nodes_points.push_back(cp_1);
					trunk_node_1.connecting_trunk_nodes_points.push_back(d11);

					if (index == 3)
					{
						
						trunk_node_0.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_exit_point, offset0));
						trunk_node_0.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_entry_point, offset0));
					}
					else
					{
						trunk_node_0.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_entry_point, offset0));
						trunk_node_0.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_exit_point, offset0));
					}


					

					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 2);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 1);

					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 2);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 1);

				}
			}


			for (int i = 0; i < trunk_node_0.connecting_trunk_nodes_points.size() && !b; i = i + 2)
			{
				double cp_0 = trunk_node_0.connecting_trunk_nodes_points[i];
				double cp_1 = trunk_node_0.connecting_trunk_nodes_points[i + 1];

				Vector2d v0 = GetOnePointFromOffset(cp_0, offset0);
				Vector2d v1 = GetOnePointFromOffset(cp_1, offset0);

				double d0 = Strip::Distance(v0, offsets_parts_0);
				double d1 = Strip::Distance(v1, offsets_parts_0);

				//bool b0 = PointsOnSharingEdge(v0, offsets_parts_0, offsets_parts_1);
				//bool b1 = PointsOnSharingEdge(v1, offsets_parts_0, offsets_parts_1);


				//if (b0&&b1)
				if (d0 < 0.00001&&d1 < 0.00001)
				{
					double d00 = DeltaDEuclideanDistance(cp_0, toolpath_size, offset0);
					Vector2d v00 = GetOnePointFromOffset(d00, offset0);

					Vector2d next_entry_point;
					Vector2d next_exit_point;

					ComputeNextEntryExitPointForInner(toolpath_size, offset0, offset1, v0, v00, next_entry_point, next_exit_point);

					b = true;

					trunk_node_0.connecting_trunk_nodes_points.push_back(d00);
					trunk_node_0.connecting_trunk_nodes_points.push_back(cp_0);

					trunk_node_1.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_exit_point, offset1));
					trunk_node_1.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_entry_point, offset1));

					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 2);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 1);

					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 2);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 1);

				}
			}


			for (int i = 0; i < trunk_node_1.connecting_trunk_nodes_points.size() && !b; i = i + 2)
			{
				double cp_0 = trunk_node_1.connecting_trunk_nodes_points[i];
				double cp_1 = trunk_node_1.connecting_trunk_nodes_points[i + 1];
				Vector2d v0 = GetOnePointFromOffset(cp_0, offset1);
				Vector2d v1 = GetOnePointFromOffset(cp_1, offset1);

				double d0 = Strip::Distance(v0, offsets_parts_1);
				double d1 = Strip::Distance(v1, offsets_parts_1);

				double smaller_d = Small_cut_direction(offset1, cp_0, cp_1);

				if (d0 < 0.00001&&d1 < 0.00001)
				{
					
					double d11 = DeltaDEuclideanDistance(smaller_d, -toolpath_size, offset1);
					Vector2d v11 = GetOnePointFromOffset(d11, offset1);

					Vector2d next_entry_point;
					Vector2d next_exit_point;

					ComputeNextEntryExitPointForOuter(toolpath_size, offset1, offset0, v1, v11, next_entry_point, next_exit_point);


					b = true;

					trunk_node_1.connecting_trunk_nodes_points.push_back(smaller_d);
					trunk_node_1.connecting_trunk_nodes_points.push_back(d11);

					trunk_node_0.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_entry_point, offset0));
					trunk_node_0.connecting_trunk_nodes_points.push_back(Circuit::FindNearestPointPar(next_exit_point, offset0));

					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 2);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_index_0);
					trunk_node_1.connecting_trunk_nodes_id.push_back(trunk_node_0.connecting_trunk_nodes_points.size() - 1);

					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 2);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_index_1);
					trunk_node_0.connecting_trunk_nodes_id.push_back(trunk_node_1.connecting_trunk_nodes_points.size() - 1);

				}
			}


			//std::vector<std::vector<Vector2d>>().swap(offset_0_pathes);
			//std::vector<std::vector<Vector2d>>().swap(offset_1_pathes);

			return b;
		}

	};

}