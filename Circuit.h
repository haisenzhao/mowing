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

		static void SelectOnePartOffset(std::vector<Vector2d> &input_points, double d, std::vector<Vector2d> &vecs)
		{
			double d1 = d + 0.00001;
			if (d1 > 1.0)
				d1 = d1 - 1.0;

			SelectOnePartOffset(input_points, d, d1, vecs);
			vecs.erase(vecs.begin() + vecs.size() - 1);
			vecs.push_back(vecs[0]);
		}

		static bool compare_double(const double &d0, const double &d1)
		{
			return d0 > d1;
		}

		static void CutACircuit(std::vector<Vector2d> &input_points, std::vector<double> &cutting_points, std::vector<std::vector<Vector2d>> &pathes)
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

					SelectOnePartOffset(input_points, cutting_points[cur_par_index + 1], cutting_points[next_par_index], one_path);

					pathes.push_back(one_path);

					std::vector<Vector2d>().swap(one_path);
				}

			}

		}

		static void SelectOnePartOffset(std::vector<Vector2d> &input_points, double d0, double d1, std::vector<Vector2d> &vecs)
		{
			if (abs(d0 - d1) < 0.0000001)
			{
				/*
				d1 = d0 + 0.00001;
				if (d1 > 1.0)
					d1 = d1 - 1.0;

				SelectOnePartOffset(input_points, d0, d1, vecs);
				vecs.erase(vecs.begin() + vecs.size() - 1);
				vecs.push_back(vecs[0]);
				*/

				Vector2d v = GetOnePointFromOffset(d0, input_points);
				vecs.push_back(v);
			}
			else
			{
				Vector2d v = GetOnePointFromOffset(d0, input_points);
				vecs.push_back(v);

				if (d1 >= 0&&d1<=1.0)
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

		static Line_2 GetRelatedLine(Vector2d v, std::vector<Vector2d> &input_points, Vector2d &v00, Vector2d &v11 )
		{
			assert(Distance(v, input_points)<0.0001);

			Line_2 line_2;

			int index = Strip::CheckSamePoint(v, input_points);

			if (index >= 0)
			{
				Vector2d v1 = input_points[(index + input_points.size() - 1) % input_points.size()];
				line_2 = Line_2(Point_2(v[0], v[1]), Point_2(v1[0], v1[1]));

				v00 = v;
				v11 = v1;
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

						v00 = v;
						v11 = v3;
						break;
					}
					length += l;
				}
			}
			return line_2;
		}
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

		static double DeltaDGeodesicDistance(double d, double distance, std::vector<Vector2d> &input_points)
		{
			double total_length = 0.0;
			for (int i = 0; i < input_points.size(); i++)
			{
				total_length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
			}

			if (distance > 0)
			{
				if (d + distance / total_length < 1.0)
				{
					return d + distance / total_length;
				}
				else
				{
					return d + distance / total_length - 1.0;
				}
			}
			else
			{
				if (d + distance / total_length > 0.0)
				{
					return d + distance / total_length;
				}
				else
				{
					return d + distance / total_length + 1.0;
				}
			}
		}


		static double Angle(Vector2d a, Vector2d b)
		{
			double ab, a1, b1, cosr;
			ab = a[0] * b[0] + a[1] * b[1];
			a1 = sqrt(a[0] * a[0] + a[1] * a[1]);
			b1 = sqrt(b[0] * b[0] + b[1] * b[1]);
			cosr = ab / a1 / b1;
			return  acos(cosr);
		}


		static void ComputeNextEntryExitPoint(double toolpath_size,std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset, Vector2d &entry_point, Vector2d &exit_point, Vector2d &next_entry_point, Vector2d &next_exit_point)
		{
			//four turning points of the first layer and second layer local_offsets
			
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
					exit_d_1 = Circuit::ComputeNextTurningPoint_complex(entry_d_1, toolpath_size, inside_offset);
				}
				else
				{
					entry_d_1 = Circuit::ComputeNextTurningPoint_complex(exit_d_1, toolpath_size, inside_offset);
				}

				entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, inside_offset);
				exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, inside_offset);
			}
			next_entry_point = entry_p_1;
			next_exit_point = exit_p_1;
		}

		static double ComputeNextTurningPoint_complex(double d, double distance, std::vector<Vector2d> &input_points)
		{
			Vector2d v = GetOnePointFromOffset(d, input_points);
			int index = Strip::CheckSamePoint(v, input_points);
			if (index >= 0)
			{
				Vector2d v0 = input_points[(index + input_points.size() - 1) % input_points.size()];
				Vector2d v1 = input_points[(index + input_points.size() + 1) % input_points.size()];

				double angle = Angle(Vector2d(v0[0] - v[0], v0[1] - v[1]), Vector2d(v1[0] - v[0], v1[1] - v[1]));

				/*
				if (angle >= PI / 2.0)
				{
					return DeltaDEuclideanDistance(d, distance, input_points);
				}
				else
				{
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
				*/

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

	};

}