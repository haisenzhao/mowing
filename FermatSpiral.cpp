#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"


namespace hpcg {


	void ToolpathGenerator::ComputeOffsets()
	{
		double lOffset = toolpath_size / 2.0;
		PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset, contours);
		while (offset_polygons.size()>0)
		{
			for (PolygonPtrVector::const_iterator pi = offset_polygons.begin(); pi != offset_polygons.end(); ++pi)
			{
				std::vector<Vector2d> one_path;
				for (Polygon_2::Vertex_const_iterator vi = (**pi).vertices_begin(); vi != (**pi).vertices_end(); ++vi)
				{
					one_path.push_back(Vector2d((*vi).x(), (*vi).y()));
				}

				offsets.push_back(one_path);
			}

			lOffset = lOffset + toolpath_size;
			offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset, contours);
		}
	}


	bool CompareTwoDouble(double d0, double d1, double d)
	{
		if (d0 < d1)
		{
			if (d >= d0 && d <= d1)
				return true;
			else
				return false;
		}
		else
		{
			if (d >= d1 && d <= d0)
				return false;
			else
				return true;
		}
	}

	void ToolpathGenerator::GenerateFermedial_axis_pointspiral()
	{
		
		ComputeOffsets();
		//BuildeImageSpace();

		bool b0 = false;

		std::vector<Vector2d> offset;
		GenerateOffset(false, offsets.size() - 1, toolpath_size / 2, offset);
		if (offset.size() == 0)
			b0=true;
		std::vector<Vector2d>().swap(offset);

		for (int i = 0; i < offsets.size(); i++)
		{
		    //double entry_d_1 = DeltaDEuclideanDistance(entry_d_0, toolpath_size, i);
			//double exit_d_1 = DeltaDEuclideanDistance(exit_d_0, toolpath_size, i);

			bool goon = true;

			double entry_d_1 = ComputeNextTurningPoint(entry_d_0, toolpath_size, i);
			double exit_d_1 = ComputeNextTurningPoint(exit_d_0, toolpath_size, i);

			bool move_entry_d_1 = true;
			bool move_exit_d_1 = true;

			if (!CompareTwoDouble(entry_d_0, exit_d_0, entry_d_1))
			{
				double d = DeltaDEuclideanDistance(entry_d_0, toolpath_size, i);

				if (!CompareTwoDouble(entry_d_0, exit_d_0, d))
				{
					entry_d_1 = exit_d_0;
				}
				else
				{
					entry_d_1 = d;
				}

				entry_d_1 = exit_d_0;
				move_entry_d_1 = false;
				//goon = false;
			}

			if (!CompareTwoDouble(exit_d_0, entry_d_0, exit_d_1))
			{
				double d = DeltaDEuclideanDistance(exit_d_0, toolpath_size, i);
				if (!CompareTwoDouble(exit_d_0, entry_d_0, d))
				{

					exit_d_1 = entry_d_0;
				}
				else
				{
					exit_d_1 = d;
				}

				exit_d_1 = entry_d_0;
				//goon = false;

				move_exit_d_1 = false;
			}

			Vector2d entry_p_0 = GetOnePointFromOffset(entry_d_0, i);
			Vector2d entry_p_1 = GetOnePointFromOffset(entry_d_1, i);
			Vector2d exit_p_0 = GetOnePointFromOffset(exit_d_0, i);
			Vector2d exit_p_1 = GetOnePointFromOffset(exit_d_1, i);


			if (i < offsets.size() - 1)
			{
				double entry_d_1_0 = FindNearestPointPar(entry_p_1, i + 1);
				double exit_d_1_0 = FindNearestPointPar(exit_p_1, i + 1);

				if (abs(entry_d_1_0 - exit_d_1_0) < 0.000001)
				{
					double t = DeltaDEuclideanDistance(exit_d_1_0, toolpath_size, i + 1);

					Vector2d v=GetOnePointFromOffset(t,i+1);
					t = FindNearestPointPar(v, i);

					if (!move_exit_d_1)
					{
						entry_d_1 = t;
					}
					else
					{
						if (!move_entry_d_1)
						{
							exit_d_1 = t;
						}
					}
				}
				exit_p_1=GetOnePointFromOffset(exit_d_1,i);
				entry_p_1=GetOnePointFromOffset(entry_d_1,i);
			}


			bool b1;
			if (i == offsets.size() - 1)
			{
				double d0, d1;
				if (entry_d_0 > exit_d_1)
				{
					d0 = entry_d_0 - exit_d_1;
				}
				else
				{
					d0 = 1.0 - (exit_d_1 - entry_d_0);
				}

				if (exit_d_0 > entry_d_1)
				{
					d1 = exit_d_0 - entry_d_1;
				}
				else
				{
					d1 = 1.0 - (entry_d_1 - exit_d_0);
				}

				if (d1 > d0)
				{
					b1 = true;
				}
				else
				{
					b1 = false;
				}
			}

			bool b2= true;

			if (i == offsets.size()-1)
			{
				Vector2d v((entry_p_1[0] + exit_p_1[0]) / 2.0, (entry_p_1[1] + exit_p_1[1]) / 2.0);

				double t = MinimalDistance(offsets[i], v);
				double t0 = FindNearestPointPar(v,i);

				if (t < 0.000001)
				{
					if (!CompareTwoDouble(entry_d_1, exit_d_1, t0))
					{
						b2 = false;
					}
				}
			}

			if (i % 2 == 0)
			{
				cccc.push_back(entry_p_0);
				dddd.push_back(entry_p_1);
				dddd.push_back(exit_p_0);
				cccc.push_back(exit_p_1);

				if (i == offsets.size() - 1)
				{
					if (b0)
					{
						entry_spiral.push_back(entry_p_0);
						exit_spiral.push_back(exit_p_0);
						break;
						if (b1)
						{
							SelectOnePartOffset(i, exit_d_0, entry_d_1, exit_spiral);
							entry_spiral.push_back(entry_p_0);
						}
						else
						{
							SelectOnePartOffset(i, entry_d_0, exit_d_1, entry_spiral);
							exit_spiral.push_back(exit_p_0);
						}
					}
					else
					{
						if (b2)
						{
							SelectOnePartOffset(i, entry_d_0, exit_d_1, entry_spiral);
							SelectOnePartOffset(i, exit_d_0, entry_d_1, exit_spiral);
						}
						else
						{
							if (b1)
							{
								SelectOnePartOffset(i, exit_d_0, entry_d_1, exit_spiral);
								entry_spiral.push_back(entry_p_0);
							}
							else
							{
								SelectOnePartOffset(i, entry_d_0, exit_d_1, entry_spiral);
								exit_spiral.push_back(exit_p_0);
							}
						}
					}

				}
				else
				{

					SelectOnePartOffset(i, entry_d_0, exit_d_1, entry_spiral);
					SelectOnePartOffset(i, exit_d_0, entry_d_1, exit_spiral);
				}

				if (i + 1 < offsets.size())
				{
					entry_d_0 = FindNearestPointPar(exit_spiral[exit_spiral.size() - 1], i + 1);
					exit_d_0 = FindNearestPointPar(entry_spiral[entry_spiral.size() - 1],i+1);
				}
			}
			else
			{
				dddd.push_back(entry_p_0);
				cccc.push_back(entry_p_1);
				cccc.push_back(exit_p_0);
				dddd.push_back(exit_p_1);

				if (i == offsets.size() - 1)
				{
					if (b0)
					{
						exit_spiral.push_back(entry_p_0);
						entry_spiral.push_back(exit_p_0);
						break;
						if (b1)
						{
							SelectOnePartOffset(i, exit_d_0, entry_d_1, entry_spiral);
							exit_spiral.push_back(entry_p_0);
						}
						else
						{
							SelectOnePartOffset(i, entry_d_0, exit_d_1, exit_spiral);
							entry_spiral.push_back(exit_p_0);
						}
					}
					else
					{
						if (b2)
						{
							SelectOnePartOffset(i, entry_d_0, exit_d_1, exit_spiral);
							SelectOnePartOffset(i, exit_d_0, entry_d_1, entry_spiral);
						}
						else
						{
							if (b1)
							{
								SelectOnePartOffset(i, exit_d_0, entry_d_1, entry_spiral);
								exit_spiral.push_back(entry_p_0);
							}
							else
							{
								SelectOnePartOffset(i, entry_d_0, exit_d_1, exit_spiral);
								entry_spiral.push_back(exit_p_0);
							}
						}
					}

				}
				else
				{
					SelectOnePartOffset(i, entry_d_0, exit_d_1, exit_spiral);
					SelectOnePartOffset(i, exit_d_0, entry_d_1, entry_spiral);
				}

				if (i + 1 < offsets.size())
				{
					entry_d_0 = FindNearestPointPar(entry_spiral[entry_spiral.size() - 1], i + 1);
					exit_d_0 = FindNearestPointPar(exit_spiral[exit_spiral.size() - 1], i + 1);
				}
			}

			if (!goon)
				break;
		}


		/*
		int i = entry_spiral.size() - 1;

		do
		{
			if (MinimalDistance(exit_spiral, entry_spiral[i]) < toolpath_size - 0.00001)
			{
				entry_spiral.erase(entry_spiral.begin() + i);
				i = entry_spiral.size() - 1;
			}
			else
			{
				break;
			}

		} while (true);


		i = exit_spiral.size() - 1;

		do
		{
			if (MinimalDistance(entry_spiral, exit_spiral[i]) < toolpath_size - 0.00001)
			{
				exit_spiral.erase(exit_spiral.begin() + i);
				i = exit_spiral.size() - 1;
			}
			else
			{
				break;
			}

		} while (true);
		*/
	}

	double ToolpathGenerator::MinimalDistanceSegments(std::vector<Vector2d> &segments, Vector2d &v)
	{
		double m_d = MAXDOUBLE;

		for (int i = 0; i < segments.size(); i = i + 2)
		{
			double d = sqrt((double)CGAL::squared_distance(Point_2(v.x, v.y), Segment_2(Point_2(segments[i].x, segments[i].y), Point_2(segments[i + 1].x, segments[i + 1].y))));
			m_d = min(m_d, d);
		}

		return m_d;
	}

	double ToolpathGenerator::MinimalDistance(std::vector<Vector2d> &vecs, Vector2d &v)
	{
		double m_d = MAXDOUBLE;

		for (int i = 0; i < vecs.size()-1; i++)
		{
			double d = sqrt((double)CGAL::squared_distance(Point_2(v.x, v.y), Segment_2(Point_2(vecs[i].x, vecs[i].y), Point_2(vecs[i+1].x, vecs[i+1].y))));

			m_d = min(m_d,d);
		}

		return m_d;

	}
	
	void ToolpathGenerator::GenerateOffset(bool direction, int offset_index, double lOffset, std::vector<Vector2d> &offset)
	{
		std::vector<Vector2d> contour;

		if (offset_index >= 0 && offset_index <= offsets.size() - 1)
		{
			for (int i = 0; i < offsets[offset_index].size(); i++)
			{
				contour.push_back(offsets[offset_index][i]);
			}
		}
		else
		{
			if (offset_index < 0)
			{
				for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
				{
					contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				}
			}
			else
			{
				assert(false);
			}
		}

		GenerateOffset(direction, contour, lOffset, offset);

		std::vector<Vector2d>().swap(contour);
	}
	void ToolpathGenerator::GenerateOffset(bool direction, std::vector<Vector2d> &contour, double lOffset, std::vector<Vector2d> &offset)
	{
		Polygon_2 polygon;

		for (int i = 0; i < contour.size(); i++)
		{
			polygon.push_back(Point_2(contour[i][0], contour[i][1]));
		}

		if (polygon.is_simple())
		{
			if (polygon.is_clockwise_oriented())
			{
				polygon.reverse_orientation();
			}

			PolygonPtrVector offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset, polygon);

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
	void ToolpathGenerator::SelectOnePartOffset(std::vector<Vector2d> &contour, double d0, double d1, std::vector<Vector2d> &vecs)
	{
		if (abs(d0 - d1) < 0.0000001)
		{
			Vector2d v=GetOnePointFromOffset(d0,contour);
			vecs.push_back(v);
		}
		else
		{
			double total_length = 0.0;
			for (int i = 0; i < contour.size(); i++)
			{
				total_length += sqrt((double)CGAL::squared_distance(Point_2(contour[i].x, contour[i].y), Point_2(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y)));
			}
			
			std::vector<double> vec_ds;
			
			vec_ds.push_back(0.0);
			double length = 0.0;
			for (int i = 0; i < contour.size(); i++)
			{
				length += sqrt((double)CGAL::squared_distance(Point_2(contour[i].x, contour[i].y), Point_2(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y)));
				vec_ds.push_back(length / total_length);
			}

			Vector2d v=GetOnePointFromOffset(d0,contour);
			vecs.push_back(v);

			if (d0 > d1)
			{
				for (int i = vec_ds.size() - 1; i >= 0; i--)
				{
					if (vec_ds[i]<d0&&vec_ds[i]>d1)
					{
						v = GetOnePointFromOffset(vec_ds[i],contour);

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
						v = GetOnePointFromOffset(vec_ds[i],contour);

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
						v = GetOnePointFromOffset(vec_ds[i],contour);
						if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
						{
							vecs.push_back(v);
						}
					}
				}
			}

			v=GetOnePointFromOffset(d1,contour);
			if (!(abs(v[0] - vecs[vecs.size() - 1][0]) < 0.000001&&abs(v[1] - vecs[vecs.size() - 1][1]) < 0.000001))
			{
				vecs.push_back(v);
			}

			std::vector<double>().swap(vec_ds);
		}

		
	}
	void ToolpathGenerator::SelectOnePartOffset(int offset_index, double d0, double d1, std::vector<Vector2d> &vecs)
	{
		std::vector<Vector2d> contour;

		if (offset_index >= 0 && offset_index <= offsets.size() - 1)
		{
			for (int i = 0; i < offsets[offset_index].size(); i++)
			{
				contour.push_back(offsets[offset_index][i]);
			}
		}
		else
		{
			if (offset_index < 0)
			{
				for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
				{
					contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				}
			}
			else
			{
				assert(false);
			}
		}

		SelectOnePartOffset(contour,d0, d1,vecs);
		std::vector<Vector2d>().swap(contour);
	}
	double ToolpathGenerator::FindNearestPointPar(Vector2d v, std::vector<Vector2d> &contour)
	{
		Vector2d n_p;

		double total_length = 0.0;
		for (int i = 0; i < contour.size(); i++)
		{
			total_length += sqrt((double)CGAL::squared_distance(Point_2(contour[i].x, contour[i].y), Point_2(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y)));
		}

		double min_d = MAXDOUBLE;
		int min_i = -1;
		for (int i = 0; i < contour.size(); i++)
		{
			Point_2 p0(contour[i].x, contour[i].y);
			Point_2 p1(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y);

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
				length += sqrt((double)CGAL::squared_distance(Point_2(contour[i].x, contour[i].y), Point_2(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y)));
			}

			Point_2 p0(contour[min_i].x, contour[min_i].y);
			Point_2 p1(contour[(min_i + 1) % contour.size()].x, contour[(min_i + 1) % contour.size()].y);

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

			length += sqrt((double)CGAL::squared_distance(Point_2(contour[min_i].x, contour[min_i].y), Point_2(n_p[0], n_p[1])));

			return length / total_length;
		}
		else
		{
			assert(false);
		}

		return -1.0;



	}
	
	double ToolpathGenerator::FindNearestPointPar(Vector2d v, int offset_index)
	{
		std::vector<Vector2d> contour;

		if (offset_index >= 0 && offset_index <= offsets.size() - 1)
		{
			for (int i = 0; i < offsets[offset_index].size(); i++)
			{
				contour.push_back(offsets[offset_index][i]);
			}
		}
		else
		{
			if (offset_index < 0)
			{
				for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
				{
					contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				}
			}
			else
			{
				assert(false);
			}
		}

		double returnD = FindNearestPointPar(v, contour);
		std::vector<Vector2d>().swap(contour);

		return returnD;
	}
	
	double ToolpathGenerator::DeltaDGeodesicDistance(double d, double distance, std::vector<Vector2d> &contour)
	{
		double total_length = 0.0;
		for (int i = 0; i < contour.size(); i++)
		{
			total_length += sqrt((double)CGAL::squared_distance(Point_2(contour[i].x, contour[i].y), Point_2(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y)));
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

	double ToolpathGenerator::ComputeNextTurningPoint(double d, double distance, int offset_index)
	{
		if (offset_index == offsets.size() - 1)
		{
			return DeltaDEuclideanDistance(d,distance,offset_index);
		}
		else
		{
	
			double d0 = DeltaDEuclideanDistance(d, distance, offset_index);
			Vector2d v0=GetOnePointFromOffset(d0,offset_index);
			double d1 = FindNearestPointPar(v0, offset_index + 1);
			Vector2d v1=GetOnePointFromOffset(d1, offset_index + 1);
			Vector2d v2=GetOnePointFromOffset(d,offset_index);
			double d2 = FindNearestPointPar(v2, offset_index + 1);
			
			if (abs(d2-d1)>0.0001)
			{
				return d0;
			}
			else
			{

				Vector2d vec(v2[0] - v1[0], v2[1] - v1[1]);
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

				Line_2 l2(Point_2(v1[0], v1[1]), Point_2(v1[0] + r_vec[0], v1[1] + r_vec[1]));

				std::vector<Vector2d> vecs;

				for (int i = 0; i < offsets[offset_index].size(); i++)
				{
					Point_2 p0(offsets[offset_index][i][0], offsets[offset_index][i][1]);
					Point_2 p1(offsets[offset_index][(i + 1) % offsets[offset_index].size()][0], offsets[offset_index][(i + 1) % offsets[offset_index].size()][1]);

					CGAL::Object result = CGAL::intersection(Segment_2(p0, p1), l2);

					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
					{
						vecs.push_back(Vector2d(ipoint->x(), ipoint->y()));
					}
				}

				for (int i = 0; i < vecs.size(); i++)
				{
					double delta_d = FindNearestPointPar(vecs[i], offset_index);

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
						std::vector<Vector2d>().swap(vecs);

						double min_d = MAXDOUBLE;
						double return_d = delta_d;
						for (int j = 1; j < 30; j++)
						{
							double t = d + j*(delta_d - d) / 30.0;

							if (t > 1.0)
								t = t - 1.0;

							Vector2d v=GetOnePointFromOffset(t,offset_index);

							double length = sqrt((double)CGAL::squared_distance(Point_2(v2[0], v2[1]), Line_2(Point_2(v[0], v[1]), Point_2(v1[0], v1[1]))));
							if (abs(length - distance) < min_d)
							{
								min_d = abs(length - distance);
								return_d = t;
							}
						}

						return return_d;
					}
				}

				assert(false);
				return -1;
			}

		}
	}

	double ToolpathGenerator::DeltaDEuclideanDistance(double d, double distance, int offset_index)
	{
		assert(d >= 0.0&&d <= 1.0);

		std::vector<Vector2d> contour;

		if (offset_index >= 0 && offset_index <= offsets.size() - 1)
		{
			for (int i = 0; i < offsets[offset_index].size(); i++)
			{
				contour.push_back(offsets[offset_index][i]);
			}
		}
		else
		{
			if (offset_index < 0)
			{
				for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
				{
					contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				}
			}
			else
			{
				assert(false);
			}
		}

		double returnd=DeltaDEuclideanDistance(d, distance, contour);
		std::vector<Vector2d>().swap(contour);
		return returnd;
	}

	double ToolpathGenerator::DeltaDEuclideanDistance(double d, double distance, std::vector<Vector2d> &contour)
	{
		Vector2d v=GetOnePointFromOffset(d,contour);

		std::vector<Vector2d> vecs;


		int divided_nb = 30;
		for (int i = 0; i < divided_nb; i++)
		{
			Point_2 p0(v[0] + abs(distance)*sin(i * 2 * PI / (double)divided_nb), v[1] + abs(distance)*cos(i * 2 * PI / (double)divided_nb));
			Point_2 p1(v[0] + abs(distance)*sin((i + 1) * 2 * PI / (double)divided_nb), v[1] + abs(distance)*cos((i + 1) * 2 * PI / (double)divided_nb));

			for (int j = 0; j < contour.size(); j++)
			{
				Point_2 p2(contour[j].x, contour[j].y);
				Point_2 p3(contour[(j + 1) % contour.size()].x, contour[(j + 1) % contour.size()].y);

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
				double delta_d = FindNearestPointPar(vecs[i], contour);

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
				double delta_d = FindNearestPointPar(vecs[i], contour);

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



	Vector2d ToolpathGenerator::GetOnePointFromOffset(double d,std::vector<Vector2d> &contour)
	{
		Vector2d v;
		double length = 0.0;
		for (int i = 0; i < contour.size(); i++)
		{
			length += sqrt((double)CGAL::squared_distance(Point_2(contour[i].x, contour[i].y), Point_2(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y)));
		}

		double total_length = length;
		length = 0.0;

		for (int i = 0; i < contour.size(); i++)
		{
			double l = sqrt((double)CGAL::squared_distance(Point_2(contour[i].x, contour[i].y), Point_2(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y)));

			if (d*total_length >= length&&d*total_length <= length + l)
			{
				double ll = (d - length / total_length)*total_length / l;
				v[0] = contour[i].x + (contour[(i + 1) % contour.size()].x - contour[i].x)*ll;
				v[1] = contour[i].y + (contour[(i + 1) % contour.size()].y - contour[i].y)*ll;
				break;
			}
			length += l;
		}

		return v;
	}

	Vector2d ToolpathGenerator::GetOnePointFromOffset(double d,int offset_index)
	{
		assert(d >= 0.0&&d <= 1.0);

		std::vector<Vector2d> contour;

		if (offset_index >= 0 && offset_index <= offsets.size() - 1)
		{
			for (int i = 0; i < offsets[offset_index].size(); i++)
			{
				contour.push_back(offsets[offset_index][i]);
			}
		}
		else
		{
			if (offset_index < 0)
			{
				for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
				{
					contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
				}
			}
			else
			{
				assert(false);
			}
		}

		Vector2d v=GetOnePointFromOffset(d,contour);
		std::vector<Vector2d>().swap(contour);
		return v;

	}

}