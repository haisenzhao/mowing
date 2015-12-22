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

		static void CuttingPoints(Vector2d &v, std::vector<Vector2d> &input_points, std::vector<Vector2d> &vecs)
		{
			std::vector<Vector2d> vt;

			double distace = Distance(v, input_points);

			for (int i = 0; i < input_points.size(); i++)
			{
				Point_2 p0(input_points[i].x, input_points[i].y);
				Point_2 p1(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y);

				double l = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), Segment_2(p0, p1)));

				if (std::abs(l - distace) < 0.00001)
				{
					double l0 = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), p0));
					double l1 = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), p1));

					if (l < l0 && l < l1)
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

						CGAL::Object result = CGAL::intersection(Line_2(p0, p1), Line_2(Point_2(v[0], v[1]), Point_2(v[0] + r_vec[0], v[1] + r_vec[1])));

						if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
						{
							bool insert = true;
							for (int j = 0; j < vt.size(); j++)
							{
								if (std::abs(vt[j][0] - ipoint->x()) < 0.000001&&std::abs(vt[j][1] - ipoint->y()) < 0.000001)
								{
									insert = false;
									break;
								}
							}
							if (insert)
							{
								vt.push_back(Vector2d(ipoint->x(), ipoint->y()));
							}
						}
						else
						{
							assert(false);
						}
					}
					else
					{
						Vector2d vp;

						if (l0 < l1)
						{
							//vt.push_back(Vector2d(p0[0], p0[1]));
							vp = Vector2d(p0[0], p0[1]);
						}
						else
						{
							//vt.push_back(Vector2d(p1[0], p1[1]));
							vp = Vector2d(p1[0], p1[1]);
						}

						bool insert = true;
						for (int j = 0; j < vt.size(); j++)
						{
							if (std::abs(vt[j][0] - vp[0]) < 0.000001&&std::abs(vt[j][1] - vp[1]) < 0.000001)
							{
								insert = false;
								break;
							}
						}
						if (insert)
						{
							vt.push_back(vp);
						}
					}
				}
			}

			for (int i = 0; i < vt.size(); i++)
			{
				vecs.push_back(vt[i]);
			}

			std::vector<Vector2d>().swap(vt);
		}

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

		static void SelectOnePartOffset(std::vector<Vector2d> &input_points, double d0, double d1, std::vector<Vector2d> &vecs)
		{
			if (abs(d0 - d1) < 0.0000001)
			{
				d1 = d0 + 0.00001;
				if (d1 > 1.0)
					d1 = d1 - 1.0;

				SelectOnePartOffset(input_points, d0, d1, vecs);
				vecs.erase(vecs.begin() + vecs.size() - 1);
				vecs.push_back(vecs[0]);

			}
			else
			{
				double total_length = 0.0;
				for (int i = 0; i < input_points.size(); i++)
				{
					total_length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
				}

				std::vector<double> vec_ds;

				vec_ds.push_back(0.0);
				double length = 0.0;
				for (int i = 0; i < input_points.size(); i++)
				{
					length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
					vec_ds.push_back(length / total_length);
				}

				Vector2d v = GetOnePointFromOffset(d0, input_points);
				vecs.push_back(v);

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

		static void DecomposeSubregions(std::vector<Vector2d> &cut_points, std::vector<Vector2d> &input_points, std::vector<std::vector<Vector2d>> &polygons)
		{
			assert(cut_points.size() % 2 == 0);

			//Decompose the whole region into sub-regions with the cutting input_points
			polygons.push_back(input_points);

			for (int i = 0; i < cut_points.size(); i = i + 2)
			{
				int polygon_index = -1;
				for (int j = 0; j < polygons.size(); j++)
				{
					double d0 = Distance(cut_points[i], polygons[j]);
					double d1 = Distance(cut_points[i + 1], polygons[j]);

					if (Distance(cut_points[i], polygons[j]) < 0.00001&&Distance(cut_points[i + 1], polygons[j]) < 0.00001)
					{
						polygon_index = j;
						break;
					}
				}



				if (polygon_index >= 0)
				{
					double d0 = FindNearestPointPar(cut_points[i], polygons[polygon_index]);
					double d1 = FindNearestPointPar(cut_points[i + 1], polygons[polygon_index]);

					std::vector<Vector2d> half_part_0;
					std::vector<Vector2d> half_part_1;

					SelectOnePartOffset(polygons[polygon_index], d0, d1, half_part_0);
					SelectOnePartOffset(polygons[polygon_index], d1, d0, half_part_1);

					polygons.erase(polygons.begin() + polygon_index);
					polygons.push_back(half_part_0);
					polygons.push_back(half_part_1);

					std::vector<Vector2d>().swap(half_part_0);
					std::vector<Vector2d>().swap(half_part_1);
				}
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

		static Vector2d ComputeEntryExitPoint(Vector2d v, double d, std::vector<Vector2d> &input_points)
		{
			std::vector<Vector2d> offset;
			GenerateOffset(input_points, d, offset);
			Vector2d vt;
			double d0 = FindNearestPointPar(v, offset);
			vt = GetOnePointFromOffset(d0, offset);
			std::vector<Vector2d>().swap(offset);
			return vt;
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

	};

}