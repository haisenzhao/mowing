#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"
#include "Strip.h"
#include "Circuit.h"

namespace hpcg {


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

	void ToolpathGenerator::FermatsSpiralTrick(std::vector<std::vector<Vector2d>> &local_offsets,
		Vector2d input_entry_point, Vector2d input_exit_point, Vector2d &output_entry_point, Vector2d &output_exit_point)
	{
		entry_d_0 = Circuit::FindNearestPointPar(input_entry_point, local_offsets[0]);
		exit_d_0 = Circuit::FindNearestPointPar(input_exit_point, local_offsets[0]);

		for (int i = 0; i < local_offsets.size() - 1; i++)
		{
			if (debug_int_1 != -100)
			{
				if (debug_int_1 >= 0)
				{
					if (i == debug_int_1)
					{
						int dsd = 0;
						break;
					}
				}

				if (debug_int_1 < 0)
				{
					if (i == abs(debug_int_1))
					{
						int dsd = 0;
					}
				}
			}


			//compute next entry and exit points
			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, local_offsets[i]);
			Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, local_offsets[i]);
			Vector2d entry_p_1, exit_p_1;
			Circuit::ComputeNextEntryExitPointForInner(toolpath_size, local_offsets[i],local_offsets[i + 1], entry_p_0, exit_p_0, entry_p_1, exit_p_1);
			double entry_d_1 = Circuit::FindNearestPointPar(entry_p_1, local_offsets[i + 1]);
			double exit_d_1 = Circuit::FindNearestPointPar(exit_p_1, local_offsets[i + 1]);

			
			std::vector<Vector2d> entry_half;
			std::vector<Vector2d> exit_half;
			Circuit::SelectOnePartOffset(local_offsets[i], entry_d_0, exit_d_0, entry_half);
			Circuit::SelectOnePartOffset(local_offsets[i], exit_d_0, entry_d_0, exit_half);

			// handle the path (entry_d_0 -> exit_d_0)
			double entry_cutting_d = Strip::IntersectPoint(entry_half, Circuit::GetRelatedLine(exit_p_1, local_offsets[i + 1]));

			if (entry_cutting_d < 0)entry_cutting_d = 0.0;
			double d0 = Circuit::DeltaDistance(exit_d_0, toolpath_size, local_offsets[i]);

			Vector2d v = Circuit::GetOnePointFromOffset(d0, local_offsets[i]);
			d0 = Strip::FindNearestPointPar(v, entry_half);
			if (entry_cutting_d>d0)
			{
				entry_cutting_d = d0;
			}
			
			// handle the path (exit_d_0 -> entry_d_0)
			double exit_cutting_d = Strip::IntersectPoint(exit_half, Circuit::GetRelatedLine(entry_p_1, local_offsets[i + 1]));
			if (exit_cutting_d < 0)exit_cutting_d = 0.0;
			d0 = Circuit::DeltaDistance(entry_d_0, toolpath_size, local_offsets[i]);
			v = Circuit::GetOnePointFromOffset(d0, local_offsets[i]);
			d0 = Strip::FindNearestPointPar(v, exit_half);
			if (exit_cutting_d>d0)
			{
				exit_cutting_d = d0;
			}

			//inpute data
			std::vector<Vector2d> path;
			Strip::SelectOnePart(entry_half, 0.0, entry_cutting_d, path);

			if (i % 2 == 0)
			{
				for (int j = 0; j < path.size(); j++)
				{
					entry_spiral.push_back(path[j]);
				}
				entry_spiral.push_back(exit_p_1);
			}
			else
			{
				for (int j = 0; j < path.size(); j++)
				{
					exit_spiral.push_back(path[j]);
				}
				exit_spiral.push_back(exit_p_1);
			}
			std::vector<Vector2d>().swap(path);

			//inpute data
			Strip::SelectOnePart(exit_half, 0.0, exit_cutting_d, path);
			if (i % 2 == 0)
			{
				for (int j = 0; j < path.size(); j++)
				{
					exit_spiral.push_back(path[j]);
				}
				exit_spiral.push_back(entry_p_1);
			}
			else
			{
				for (int j = 0; j < path.size(); j++)
				{
					entry_spiral.push_back(path[j]);
				}
				entry_spiral.push_back(entry_p_1);
			}

			std::vector<Vector2d>().swap(path);
			std::vector<Vector2d>().swap(entry_half);
			std::vector<Vector2d>().swap(exit_half);

			entry_d_0 = entry_d_1;
			exit_d_0 = exit_d_1;
		}

		output_entry_point = Circuit::GetOnePointFromOffset(entry_d_0, local_offsets[local_offsets.size() - 1]);
		output_exit_point = Circuit::GetOnePointFromOffset(exit_d_0, local_offsets[local_offsets.size() - 1]);

		if (debug_int_1 != -100)
		{
			return;
		}

		//final offset
		#pragma region final_offset
		Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, local_offsets[local_offsets.size() - 1]);
		Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, local_offsets[local_offsets.size() - 1]);

		double entry_d_1 = Circuit::DeltaDistance(exit_d_0, toolpath_size, local_offsets[local_offsets.size() - 1]);
		double exit_d_1 = Circuit::DeltaDistance(entry_d_0, toolpath_size, local_offsets[local_offsets.size() - 1]);

		Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, local_offsets[local_offsets.size() - 1]);
		Vector2d exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, local_offsets[local_offsets.size() - 1]);


		std::vector<Vector2d> vecs0, vecs1;
		Circuit::SelectOnePartOffset(local_offsets[local_offsets.size() - 1], entry_d_0, entry_d_1, vecs0);
		Circuit::SelectOnePartOffset(local_offsets[local_offsets.size() - 1], exit_d_0, exit_d_1, vecs1);

	
		if ((local_offsets.size() - 1) % 2 == 0)
		{
			if (Strip::GetTotalLength(vecs0) > Strip::GetTotalLength(vecs1))
			{
				for (int i = 1; i < vecs0.size(); i++)
				{
					entry_spiral.push_back(vecs0[i]);
				}
			}
			else
			{
				for (int i = 1; i < vecs1.size(); i++)
				{
					exit_spiral.push_back(vecs1[i]);
				}

			}
		}
		else
		{
			if (Strip::GetTotalLength(vecs0) > Strip::GetTotalLength(vecs1))
			{
				for (int i = 1; i < vecs0.size(); i++)
				{
					exit_spiral.push_back(vecs0[i]);
				}
			}
			else
			{
				for (int i = 1; i < vecs1.size(); i++)
				{
					entry_spiral.push_back(vecs1[i]);
				}
			}
		}



		std::vector<Vector2d>().swap(vecs0);
		std::vector<Vector2d>().swap(vecs1);
		#pragma endregion

	}

	void ToolpathGenerator::FermatSpiral(std::vector<Vector2d> &contour, Vector2d input_entry_point, Vector2d input_exit_point)
	{
		Circuit::ComputeOffsets(toolpath_size,contour,offsets);

		entry_d_0 = Circuit::FindNearestPointPar(input_entry_point, offsets[0]);
		exit_d_0 = Circuit::FindNearestPointPar(input_exit_point, offsets[0]);

		for (int i = 0; i < offsets.size(); i++)
		{
			
			
			double entry_d_1 = Circuit::DeltaDEuclideanDistance(entry_d_0, toolpath_size, offsets[i]);
			double exit_d_1 = Circuit::DeltaDEuclideanDistance(exit_d_0, toolpath_size, offsets[i]);
			
			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offsets[i]);
			Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, offsets[i]);
			Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, offsets[i]);
			Vector2d exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, offsets[i]);

			//check vilid of entry and exit points
			bool change_entry_d_1 = false;
			bool change_exit_d_1 = false;
			if (!CompareTwoDouble(entry_d_0, exit_d_0, entry_d_1))
			{
				entry_d_1 = exit_d_0;
				change_entry_d_1 = true;
			}

			if (!CompareTwoDouble(exit_d_0, entry_d_0, exit_d_1))
			{
				exit_d_1 = entry_d_0;
				change_exit_d_1 = true;
			}

			//handle the case where the next_exit_d_0 and next_entry_d_0 meets each other
			if (i < offsets.size() - 1)
			{
				double next_exit_d_0 = Circuit::FindNearestPointPar(entry_p_1, offsets[i + 1]);
				double next_entry_d_0 = Circuit::FindNearestPointPar(exit_p_1, offsets[i + 1]);

				if (abs(next_exit_d_0 - next_entry_d_0) < 0.000001)
				{
					if (i % 2 == 0)
					{
						if (!change_exit_d_1)
						{
							double next_entry_d_1 = Circuit::DeltaDEuclideanDistance(next_entry_d_0, toolpath_size, offsets[i + 1]);
							Vector2d next_entry_p_1 = Circuit::GetOnePointFromOffset(next_entry_d_1, offsets[i + 1]);
							exit_d_1 = Circuit::FindNearestPointPar(next_entry_p_1, offsets[i]);
						}
						else
						{
							if (!change_entry_d_1)
							{
								double next_exit_d_1 = Circuit::DeltaDEuclideanDistance(next_exit_d_0, toolpath_size, offsets[i + 1]);
								Vector2d next_entry_p_1 = Circuit::GetOnePointFromOffset(next_exit_d_1, offsets[i + 1]);
								entry_d_1 = Circuit::FindNearestPointPar(next_entry_p_1, offsets[i]);
							}
						}
					}
					else
					{
						if (!change_entry_d_1)
						{
							double next_exit_d_1 = Circuit::DeltaDEuclideanDistance(next_exit_d_0, toolpath_size, offsets[i + 1]);
							Vector2d next_entry_p_1 = Circuit::GetOnePointFromOffset(next_exit_d_1, offsets[i + 1]);
							entry_d_1 = Circuit::FindNearestPointPar(next_entry_p_1, offsets[i]);
						}
						else
						{
							if (!change_exit_d_1)
							{
								double next_entry_d_1 = Circuit::DeltaDEuclideanDistance(next_entry_d_0, toolpath_size, offsets[i + 1]);
								Vector2d next_entry_p_1 = Circuit::GetOnePointFromOffset(next_entry_d_1, offsets[i + 1]);
								exit_d_1 = Circuit::FindNearestPointPar(next_entry_p_1, offsets[i]);
							}
						}
					}
				}
			}

			entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offsets[i]);
			entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, offsets[i]);
			exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, offsets[i]);
			exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, offsets[i]);

			turning_points_entry.push_back(entry_p_0);
			turning_points_exit.push_back(entry_p_1);
			turning_points_exit.push_back(exit_p_0);
			turning_points_entry.push_back(exit_p_1);

			//select part
			if (i % 2 == 0)
			{
				Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, entry_spiral);
				Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, exit_spiral);

				if (i + 1 < offsets.size())
				{
					entry_d_0 = Circuit::FindNearestPointPar(exit_spiral[exit_spiral.size() - 1], offsets[i + 1]);
					exit_d_0 = Circuit::FindNearestPointPar(entry_spiral[entry_spiral.size() - 1], offsets[i + 1]);
				}
			}
			else
			{
				Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, exit_spiral);
				Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, entry_spiral);

				if (i + 1 < offsets.size())
				{
					entry_d_0 = Circuit::FindNearestPointPar(entry_spiral[entry_spiral.size() - 1], offsets[i + 1]);
					exit_d_0 = Circuit::FindNearestPointPar(exit_spiral[exit_spiral.size() - 1], offsets[i + 1]);
				}
			}
		}




	}

	double ToolpathGenerator::ComputeNextTurningPoint(double d, double distance, int offset_index)
	{
		if (offset_index == offsets.size() - 1)
		{
			return Circuit::DeltaDEuclideanDistance(d,distance,offsets[offset_index]);
		}
		else
		{
			double d0 = Circuit::DeltaDEuclideanDistance(d, distance, offsets[offset_index]);
			Vector2d v0 = Circuit::GetOnePointFromOffset(d0, offsets[offset_index]);
			double d1 = Circuit::FindNearestPointPar(v0, offsets[offset_index + 1]);
			Vector2d v1 = Circuit::GetOnePointFromOffset(d1, offsets[offset_index + 1]);
			Vector2d v2 = Circuit::GetOnePointFromOffset(d, offsets[offset_index]);
			double d2 = Circuit::FindNearestPointPar(v2, offsets[offset_index + 1]);
			
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
					double delta_d = Circuit::FindNearestPointPar(vecs[i], offsets[offset_index]);

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

							Vector2d v = Circuit::GetOnePointFromOffset(t, offsets[offset_index]);

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

}