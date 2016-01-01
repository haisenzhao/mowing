#pragma once

#include <iostream>
#include <fstream>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>

#include <ToolpathGenerator.h>

namespace hpcg {

	class Strip
	{
	public:
		static Vector2d GetOnePointFromStrip(double d, std::vector<Vector2d> &input_points)
		{
			Vector2d v;

			double length = 0.0;
			for (int i = 0; i < input_points.size() - 1; i++)
			{
				length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
			}
			double total_length = length;

			length = 0.0;

			for (int i = 0; i < input_points.size() - 1; i++)
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

		static double FindNearestPointPar(Vector2d v, std::vector<Vector2d> &input_points)
		{
			Vector2d n_p;

			double total_length = GetTotalLength(input_points);

			double min_d = MAXDOUBLE;
			int min_i = -1;
			for (int i = 0; i < input_points.size() - 1; i++)
			{
				Point_2 p0(input_points[i].x, input_points[i].y);
				Point_2 p1(input_points[i + 1].x, input_points[i + 1].y);

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
					length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[i + 1].x, input_points[i + 1].y)));
				}

				Point_2 p0(input_points[min_i].x, input_points[min_i].y);
				Point_2 p1(input_points[min_i + 1].x, input_points[min_i + 1].y);

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

		static double Distance(Vector2d v0, Vector2d v1)
		{
			return std::sqrt((v0[0] - v1[0])*(v0[0] - v1[0]) + (v0[1] - v1[1])*(v0[1] - v1[1]));
		}

		static double Distance(Vector2d v, std::vector<Vector2d> &input_points)
		{
			double m_d = MAXDOUBLE;
			for (int i = 0; i < input_points.size()-1; i++)
			{
				double d = sqrt((double)CGAL::squared_distance(Point_2(v.x, v.y), Segment_2(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y))));
				m_d = min(m_d, d);
			}
			return m_d;
		}


		static double GetTotalLength(std::vector<Vector2d> &input_points)
		{
			double length = 0.0;

			if (input_points.size() >= 2)
			{
				for (int i = 0; i < input_points.size() - 1; i++)
				{
					length += sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(input_points[(i + 1) % input_points.size()].x, input_points[(i + 1) % input_points.size()].y)));
				}
			}

			return length;
		}

		static int  CheckSamePoint(Vector2d v, std::vector<Vector2d> &input_points)
		{
			int index = -1;
			for (int i = 0; i < input_points.size() && index<0; i++)
			{
				double l = sqrt((double)CGAL::squared_distance(Point_2(input_points[i].x, input_points[i].y), Point_2(v[0], v[1])));

				if (l < 0.000001)
				{
					index = i;
				}
			}
			return index;
		}

		static void SamplingPoints(std::vector<Vector2d> &input_points, double delta, Vector2d v, std::vector<Vector2d> &output_points)
		{
			double insert_par = FindNearestPointPar(v, input_points);
			bool b = true;
			for (double d = 0.0; d <= 1.0; d = d + delta)
			{
				output_points.push_back(GetOnePointFromStrip(d, input_points));

				if (d <= insert_par&&insert_par <= d + delta&&b)
				{
					output_points.push_back(v);
					b = false;
				}
			}
		}

		static void StripDeformation(std::vector<Vector2d> &input_points_0, std::vector<Vector2d> &input_points_1)
		{
			for (int i = 0; i < input_points_0.size(); i++)
			{
				//double t = FindNearestPointPar(input_points_0[i], input_points_0);
				double t = (double)i / (double)(input_points_0.size() - 1);
				Vector2d n(input_points_1[i][0] - input_points_0[i][0], input_points_1[i][1] - input_points_0[i][1]);
			
				input_points_0[i][0] = input_points_0[i][0] + n[0] * t;
				input_points_0[i][1] = input_points_0[i][1] + n[1] * t;
			}
		}

		static void SelectOnePart(std::vector<Vector2d> &input_points, double d0, double d1, std::vector<Vector2d> &vecs)
		{
			if (abs(d0 - d1) > 0.00001&&d1 > d0)
			{
				vecs.push_back(GetOnePointFromStrip(d0, input_points));

				for (int i = 1; i < input_points.size() - 1; i++)
				{
					double par = FindNearestPointPar(input_points[i], input_points);
					if (par > d0&&par < d1)
					{
						vecs.push_back(input_points[i]);
					}
				}
				Vector2d v = GetOnePointFromStrip(d1, input_points);

				if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
				{
					vecs.push_back(v);
				}

				if (abs(vecs[1][0] - vecs[0][0]) < 0.000001&&abs(vecs[1][1] - vecs[0][1]) < 0.000001)
				{
					vecs.erase(vecs.begin());
				}
			}

			if (abs(d0 - d1) < 0.00001)
			{
				vecs.push_back(GetOnePointFromStrip(d0, input_points));
			}
		}

		static bool TwoSegmentIntersect(Vector2d v0, Vector2d v1, Vector2d v2, Vector2d v3)
		{
			CGAL::Object result = CGAL::intersection(Segment_2(Point_2(v0[0], v0[1]), Point_2(v1[0], v1[1])), Segment_2(Point_2(v2[0], v2[1]), Point_2(v3[0], v3[1])));

			if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		static double IntersectPoint(std::vector<Vector2d> &input_points, Vector2d v0, Vector2d v1)
		{
			bool b = false;
			Vector2d v;
			//for (int i = 0; i < input_points.size() - 1&&!b; i++)
			for (int i = input_points.size() - 2; i >=0 && !b; i--)
			{
				CGAL::Object result = CGAL::intersection(Segment_2(Point_2(input_points[i][0], input_points[i][1]), Point_2(input_points[i + 1][0], input_points[i + 1][1])),
					Line_2(Point_2(v0[0], v0[1]), Point_2(v1[0], v1[1])));

				if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
				{
					v[0] = ipoint->x();
					v[1] = ipoint->y();
					b = true;
				}
			}

			if (b)
			{
				return FindNearestPointPar(v,input_points);
			}
			else
			{
				return -1.0;
			}
		}

		static double IntersectPoint(std::vector<Vector2d> &input_points, Line_2 &line_2)
		{
			bool b = false;
			Vector2d v;
			for (int i = input_points.size() - 2; i >= 0 && !b; i--)
			{
				CGAL::Object result = CGAL::intersection(Segment_2(Point_2(input_points[i][0], input_points[i][1]), Point_2(input_points[i + 1][0], input_points[i + 1][1])),line_2);

				if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
				{
					v[0] = ipoint->x();
					v[1] = ipoint->y();
					b = true;
				}
			}

			if (b)
			{
				return FindNearestPointPar(v, input_points);
			}
			else
			{
				return -1.0;
			}
		}

	};


}
