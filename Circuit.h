#pragma once

#include <iostream>
#include <fstream>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>

#include <ToolpathGenerator.h>
#include <Strip.h>

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

			if (a[0] * b[1] - a[1] * b[0] > 0)
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
			double m_d = MAXDOUBLE;
			for (int i = 0; i < input_points.size(); i++)
			{
				double d = sqrt((double)CGAL::squared_distance(Point_2(v.x, v.y), Segment_2(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y))));
				m_d = min(m_d, d);
			}
			return m_d;
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

		static void ComputeOffsets(double toolpath_size, Polygon_with_holes &contours, 
			std::vector<std::vector<Vector2d>> &offsets, std::vector<std::vector<std::vector<Vector2d>>> &offsetses)
		{
			double lOffset = toolpath_size / 2.0;
			std::vector<std::vector<Vector2d>> one_pathes;
			Circuit::GenerateOffsetHole(contours, lOffset, one_pathes);

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
				Circuit::GenerateOffsetHole(contours, lOffset, one_pathes);
			}



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

				if (min_d<l0 &&min_d<l1)
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

				Vector2d v1,v2;
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
			}
			else
			{
				return DeltaDEuclideanDistance(d, distance, input_points);
			}
		}

		static void FindOptimalEntryExitPoints(double toolpath_size,std::vector<Vector2d> &offset, Vector2d &entry_point, Vector2d &exit_point)
		{
			double m_d = MAXDOUBLE;
			int m_d_index = -1;

			for (int i = 0; i < offset.size(); i++)
			{
				Vector2d v0 = offset[(i - 1 + offset.size()) % offset.size()];
				Vector2d v1 = offset[i];
				Vector2d v2 = offset[(i + 1 + offset.size()) % offset.size()];

				double angle = Circuit::Angle2PI(Vector2d(v0[0] - v1[0], v0[1] - v1[1]), Vector2d(v2[0] - v1[0], v2[1] - v1[1]));
				if (angle <= PI)
				{
					if (abs(angle - PI / 2.0) < m_d)
					{
						m_d = abs(angle - PI / 2.0);
						m_d_index = i;
					}
				}
			}

			entry_point = offset[m_d_index];
			double d = Circuit::FindNearestPointPar(entry_point, offset);
			double entry_point_d = Circuit::DeltaDistance(d, toolpath_size, offset);
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
				double d = Strip::Distance(v0,vecs[i]);

				if (d < min_d)
				{
					min_d = d;
					index = i;
				}
			}

			return vecs[index];
		}

		static void ComputeNextEntryExitPointForOuter(double toolpath_size, std::vector<Vector2d> &outside_offset,std::vector<Vector2d> &inside_offset, 
			Vector2d &entry_point, Vector2d &exit_point, Vector2d &next_entry_point, Vector2d &next_exit_point)
		{

			Vector2d n(exit_point[0] - entry_point[0], exit_point[1] - entry_point[1]);
			Vector2d pn(-n[1],n[0]);

			next_entry_point = DetectNearestPoint(entry_point, Vector2d(entry_point[0] + pn[0], entry_point[1] + pn[1]), inside_offset);
			next_exit_point = DetectNearestPoint(exit_point, Vector2d(exit_point[0] + pn[0], exit_point[1] + pn[1]), inside_offset);
		}



		static void ComputeNextEntryExitPointForInner(double toolpath_size, std::vector<Vector2d> &outside_offset,std::vector<Vector2d> &inside_offset,
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
					double temp = Circuit::DeltaDistance(entry_d_1, toolpath_size, inside_offset);

					if (temp >= 0)
						exit_d_1 = temp;
				}
				else
				{
					double temp = Circuit::DeltaDistance(exit_d_1, toolpath_size, inside_offset);
					if (temp >= 0)
						entry_d_1 = temp;
				}

				entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, inside_offset);
				exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, inside_offset);
			}
			next_entry_point = entry_p_1;
			next_exit_point = exit_p_1;
		}


		//for fermat spiral generation
		static Line_2 GetRelatedLine(Vector2d v, std::vector<Vector2d> &input_points, std::vector<Vector2d> &debug_points)
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

			if (index>=0)
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
				Vector2d v1 = Circuit::GetOnePointFromOffset(trunk_node.connecting_points[i+1], input_points);
				Vector2d v3((v0[0] + v1[0]) / 2.0, (v0[1] + v1[1]) / 2.0);
				double par = Circuit::FindNearestPointPar(v3, input_points);

				cutting_points_par.push_back(par);
			}

			std::sort(cutting_points_par.begin(), cutting_points_par.end(), compare_double);


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
	

		static bool SharingPart(std::vector<Vector2d> &main_offset, std::vector<Vector2d> &check_offset, double &offset_start, double &offset_end, std::vector<std::vector<Vector2d>> &pathes_temp)
		{
			int index = -1;

			for (int i = 0; i < main_offset.size(); i++)
			{
				double d0 = Distance(main_offset[i], check_offset);
				double d1 = Distance(main_offset[(i + 1) % main_offset.size()], check_offset);

				if (d0 > 0.00001&&d1<0.00001)
				{
					index = (i + 1) % main_offset.size();
					break;
				}
			}



			if (index < 0)
			{
				offset_start = 0.0;
				offset_end = 1.0;
				return true;
			}

			offset_start = FindNearestPointPar(main_offset[index], main_offset);

			int i = (index + 1) % main_offset.size();
			do
			{
				double d = Distance(main_offset[i], check_offset);

				if (d > 0.00001||i==index)
				{
					break;
				}
				i = (i + 1) % main_offset.size();
			} while (true);
			
			i = (i + main_offset.size()- 1) % main_offset.size();

			offset_end = FindNearestPointPar(main_offset[i], main_offset);

			return false;
		}

		static void DetectSharingPartEnclose(double toolpath_size, std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset,
			double &outside_offset_start, double &outside_offset_end, double &inside_offset_start, double &inside_offset_end, std::vector<std::vector<Vector2d>> &pathes_temp)
		{
			std::vector<Vector2d> offset;
			GenerateOffsetForOutside(inside_offset, toolpath_size, offset);

			bool b = SharingPart(outside_offset, offset, outside_offset_start, outside_offset_end, pathes_temp);

			if (b)
			{
				inside_offset_start = 0.0;
				inside_offset_end = 1.0;
			}
			else
			{
				Vector2d v0 = GetOnePointFromOffset(outside_offset_start, outside_offset);
				Vector2d v1 = GetOnePointFromOffset(outside_offset_end, outside_offset);

				inside_offset_start = FindNearestPointPar(v0,inside_offset);
				inside_offset_end = FindNearestPointPar(v1, inside_offset);
			}

			std::vector<Vector2d>().swap(offset);
		}

		static void DetectSharingPart(double toolpath_size, std::vector<Vector2d> &offset_0, std::vector<Vector2d> &offset_1,
			double &offset_0_start, double &offset_0_end, double &offset_1_start, double &offset_1_end, std::vector<std::vector<Vector2d>> &pathes_temp)
		{
			bool b0 = CheckInside(offset_0[0], offset_1);
			bool b1 = CheckInside(offset_1[0], offset_0);

			if (b0)
			{
				DetectSharingPartEnclose(toolpath_size, offset_1, offset_0, offset_1_start, offset_1_end, offset_0_start, offset_0_end, pathes_temp);
			}

			if (b1)
			{
				DetectSharingPartEnclose(toolpath_size, offset_0, offset_1, offset_0_start, offset_0_end, offset_1_start, offset_1_end, pathes_temp);
			}

			if (!b0&&!b1)
			{
				std::vector<Vector2d> offset;
				GenerateOffsetForOutside(offset_0, toolpath_size, offset);

				bool b = SharingPart(offset_1, offset, offset_1_start, offset_1_end, pathes_temp);

				Vector2d v0 = GetOnePointFromOffset(offset_1_start, offset_1);
				Vector2d v1 = GetOnePointFromOffset(offset_1_end, offset_1);

				offset_0_start = FindNearestPointPar(v0, offset_0);
				offset_0_end = FindNearestPointPar(v1, offset_0);

				std::vector<Vector2d>().swap(offset);
			}
		}

		static void ConnectTwoTrunkNodes(double toolpath_size, std::vector<Vector2d> &offset0, std::vector<Vector2d> &offset1,
			int trunk_node_index_0, int trunk_node_index_1,
			TrunkNode &trunk_node_0, TrunkNode &trunk_node_1,
			std::vector<std::vector<Vector2d>> &pathes, std::vector<std::vector<Vector2d>> &pathes_temp, std::vector<Vector2d> &debug_points)
		{
			double offset_0_start;
			double offset_0_end;
			double offset_1_start;
			double offset_1_end;

			DetectSharingPart(toolpath_size, offset0, offset1, offset_0_start, offset_0_end, offset_1_start, offset_1_end, pathes_temp);

			std::vector<Vector2d> offset_0_path;
			SelectOnePartOffset(offset0, offset_0_end, offset_0_start, offset_0_path);
			//pathes_temp.push_back(offset_0_path);
	
			std::vector<Vector2d> offset_1_path;
			SelectOnePartOffset(offset1, offset_1_end, offset_1_start, offset_1_path);
			//pathes_temp.push_back(offset_1_path);

			bool b = false;

			for (int i = 0; i < trunk_node_0.connecting_leaf_nodes_points.size()&&!b; i = i + 2)
			{
				double cp_0 = trunk_node_0.connecting_leaf_nodes_points[i];
				double cp_1 = trunk_node_0.connecting_leaf_nodes_points[i + 1];

				Vector2d v0 = GetOnePointFromOffset(cp_0, offset0);
				Vector2d v1 = GetOnePointFromOffset(cp_1, offset0);

				double d0 = Strip::Distance(v0, offset_0_path);
				double d1 = Strip::Distance(v1, offset_0_path);

				if (d0 < 0.00001&&d1 < 0.00001)
				{
					double d00 = DeltaDEuclideanDistance(cp_0, toolpath_size, offset0);
					Vector2d v00 = GetOnePointFromOffset(d00, offset0);

					Vector2d next_entry_point;
					Vector2d next_exit_point;

					ComputeNextEntryExitPointForInner(toolpath_size, offset0, offset1, v0, v00, next_entry_point, next_exit_point);
					
					//turning_points_entry_temp.push_back(v0);
					//turning_points_exit_temp.push_back(v00);

					//turning_points_entry_temp.push_back(next_entry_point);
					//turning_points_exit_temp.push_back(next_exit_point);
					b = true;

					//cutting_points_0.push_back(d00);
					//cutting_points_0.push_back(cp_0);

					//cutting_points_1.push_back(Circuit::FindNearestPointPar(next_exit_point, offset1));
					//cutting_points_1.push_back(Circuit::FindNearestPointPar(next_entry_point, offset1));
					
					/*
					std::vector<Vector2d> one_path;
					one_path.push_back(v0);
					one_path.push_back(next_entry_point);
					pathes.push_back(one_path);
					std::vector<Vector2d>().swap(one_path);
					one_path.push_back(v00);
					one_path.push_back(next_exit_point);
					pathes.push_back(one_path);
					*/

					//int cutting_index_0, int cutting_index_1,
					//std::vector<int> &cutting_points_index_0, std::vector<int> &cutting_points_index_1,

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
				double d0 = Strip::Distance(v0, offset_1_path);
				double d1 = Strip::Distance(v1, offset_1_path);
				
				if (d0 < 0.00001&&d1 < 0.00001)
				{

					double d11 = DeltaDEuclideanDistance(cp_1, -toolpath_size, offset1);
					Vector2d v11 = GetOnePointFromOffset(d11, offset1);

					//turning_points_entry_temp.push_back(v1);
					//turning_points_exit_temp.push_back(v11);

					Vector2d next_entry_point;
					Vector2d next_exit_point;

					ComputeNextEntryExitPointForOuter(toolpath_size, offset1, offset0, v1, v11, next_entry_point, next_exit_point);
					//turning_points_entry_temp.push_back(next_entry_point);
					//turning_points_exit_temp.push_back(next_exit_point);

					b = true;


					//debug_points.push_back(v1);
					//debug_points.push_back(v11);
					//debug_points.push_back(next_entry_point);
					//debug_points.push_back(next_exit_point);
					
					
					//cutting_points_1.push_back(cp_1);
					//cutting_points_1.push_back(d11);

					//cutting_points_0.push_back(Circuit::FindNearestPointPar(next_entry_point, offset0));
					//cutting_points_0.push_back(Circuit::FindNearestPointPar(next_exit_point, offset0));

					/*
					std::vector<Vector2d> one_path;
					one_path.push_back(v1);
					one_path.push_back(next_entry_point);
					pathes.push_back(one_path);
					std::vector<Vector2d>().swap(one_path);
					one_path.push_back(v11);
					one_path.push_back(next_exit_point);
					pathes.push_back(one_path);
					*/

					//int cutting_index_0, int cutting_index_1,
					//std::vector<int> &cutting_points_index_0, std::vector<int> &cutting_points_index_1,



					trunk_node_1.connecting_trunk_nodes_points.push_back(cp_1);
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

			std::vector<Vector2d>().swap(offset_0_path);
			std::vector<Vector2d>().swap(offset_1_path);
		}

	};

}