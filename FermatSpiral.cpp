#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"
#include "Strip.h"
#include "Circuit.h"

namespace hpcg {

	void ToolpathGenerator::ComputeOffsets(std::vector<Vector2d> &contour)
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

	void ToolpathGenerator::FermatsSpiralTrick(std::vector<Vector2d> &contour, Vector2d input_entry_point, Vector2d input_exit_point)
	{
		ComputeOffsets(contour);

		entry_d_0 = Circuit::FindNearestPointPar(input_entry_point, offsets[0]);
		exit_d_0 = Circuit::FindNearestPointPar(input_exit_point, offsets[0]);

		for (int i = 0; i < offsets.size() - 1; i++)
		{
			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offsets[i]);
			Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, offsets[i]);

			double entry_d_1 = Circuit::FindNearestPointPar(entry_p_0, offsets[i + 1]);
			double exit_d_1 = Circuit::FindNearestPointPar(exit_p_0, offsets[i + 1]);

			Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, offsets[i + 1]);
			Vector2d exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, offsets[i+1]);

			//////
			std::vector<Vector2d> half;
			std::vector<Vector2d> next_half;
			Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_0, half);
			Circuit::SelectOnePartOffset(offsets[i + 1], exit_d_1, entry_d_1, next_half);
			
			//double d = Strip::IntersectPoint(half, Circuit::GetTangent(exit_p_1, offsets[i + 1]));
			double d = Strip::IntersectPoint(half, next_half[0], next_half[1]);
			double d0 = ComputeNextTurningPoint(exit_d_0, toolpath_size, i);

			if (Strip::Distance(exit_p_0, Circuit::GetOnePointFromOffset(d, offsets[i])) <
				Strip::Distance(exit_p_0, Circuit::GetOnePointFromOffset(d0, offsets[i])))
			{
				d = d0;
			}

			std::vector<Vector2d> path;
			Strip::SelectOnePart(half, 0.0, d, path);

			if (i % 2 == 0)
			{
				for (int j = 0; j < path.size(); j++)
				{
					entry_spiral.push_back(path[j]);
				}
				entry_spiral.push_back(next_half[0]);
			}
			else
			{
				for (int j = 0; j < path.size(); j++)
				{
					exit_spiral.push_back(path[j]);
				}
				exit_spiral.push_back(next_half[0]);
			}

			std::vector<Vector2d>().swap(path);
			std::vector<Vector2d>().swap(half);
			std::vector<Vector2d>().swap(next_half);

			//////
			Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_0, half);
			Circuit::SelectOnePartOffset(offsets[i + 1], entry_d_1, exit_d_1, next_half);

			//d = Strip::IntersectPoint(half, Circuit::GetTangent(entry_p_1, offsets[i + 1]));
			d = Strip::IntersectPoint(half, next_half[0], next_half[1]);

			d0 = ComputeNextTurningPoint(entry_d_0, toolpath_size, i);

			if (Strip::Distance(entry_p_0, Circuit::GetOnePointFromOffset(d, offsets[i])) <
				Strip::Distance(entry_p_0, Circuit::GetOnePointFromOffset(d0, offsets[i])))
			{
				d = d0;
			}

			Strip::SelectOnePart(half, 0.0, d, path);

			if (i % 2 == 0)
			{
				for (int j = 0; j < path.size(); j++)
				{
					exit_spiral.push_back(path[j]);
				}
				exit_spiral.push_back(next_half[0]);
			}
			else
			{
				for (int j = 0; j < path.size(); j++)
				{
					entry_spiral.push_back(path[j]);
				}
				entry_spiral.push_back(next_half[0]);
			}

			std::vector<Vector2d>().swap(path);
			std::vector<Vector2d>().swap(half);
			std::vector<Vector2d>().swap(next_half);

			entry_d_0 = entry_d_1;
			exit_d_0 = exit_d_1;

		}
	}

	void ToolpathGenerator::FermatsSpiralSmooth1(std::vector<Vector2d> &contour, Vector2d input_entry_point, Vector2d input_exit_point)
	{
		ComputeOffsets(contour);

		std::vector<Vector2d> offset_0;
		std::vector<Vector2d> offset_1;

		Circuit::GenerateOffset(contour, toolpath_size*0.5, offset_0);
		Circuit::GenerateOffset(contour, toolpath_size*1.5, offset_1);

		entry_d_0 = Circuit::FindNearestPointPar(input_entry_point, offset_0);
		exit_d_0 = Circuit::FindNearestPointPar(input_exit_point, offset_0);


		//first
		Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offset_0);
		Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, offset_0);

		double entry_d_1 = Circuit::FindNearestPointPar(entry_p_0, offset_1);
		double exit_d_1 = Circuit::FindNearestPointPar(exit_p_0, offset_1);

		//////
		std::vector<Vector2d> half;
		std::vector<Vector2d> next_half;
		Circuit::SelectOnePartOffset(offset_0, entry_d_0, exit_d_0, half);
		Circuit::SelectOnePartOffset(offset_1, exit_d_1, entry_d_1, next_half);

		std::vector<Vector2d> entry_path;
		Strip::SelectOnePart(half, 0.0, Strip::IntersectPoint(half, next_half[0], next_half[1]), entry_path);
		entry_path.push_back(next_half[0]);

		for (int j = 0; j < entry_path.size(); j++)
		{
			entry_spiral.push_back(entry_path[j]);
		}

		//std::vector<Vector2d>().swap(entry_path);
		std::vector<Vector2d>().swap(half);
		std::vector<Vector2d>().swap(next_half);

		//////
		Circuit::SelectOnePartOffset(offset_0, exit_d_0, entry_d_0, half);
		Circuit::SelectOnePartOffset(offset_1, entry_d_1, exit_d_1, next_half);

		std::vector<Vector2d> exit_path;
		Strip::SelectOnePart(half, 0.0, Strip::IntersectPoint(half, next_half[0], next_half[1]), exit_path);
		exit_path.push_back(next_half[0]);
		for (int j = 0; j < exit_path.size(); j++)
		{
			exit_spiral.push_back(exit_path[j]);
		}

		//std::vector<Vector2d>().swap(exit_path);
		std::vector<Vector2d>().swap(half);
		std::vector<Vector2d>().swap(next_half);


		//entry_path->halp
		Circuit::SelectOnePartOffset(offset_0, exit_d_0, entry_d_0, half);
		Circuit::SelectOnePartOffset(offset_0, entry_d_0, exit_d_0, next_half);
		std::vector<Vector2d>().swap(offset_0);
		std::vector<Vector2d>().swap(offset_1);

		for (int i = 1; i < entry_path.size(); i++)
		{
			half.push_back(entry_path[i]);
		}
		Circuit::GenerateOffset(half, toolpath_size, offset_0);
	
		//exit_path->next_half
		for (int i = 1; i < exit_path.size(); i++)
		{
			next_half.push_back(exit_path[i]);
		}
		Circuit::GenerateOffset(next_half, toolpath_size, offset_1);

		std::vector<Vector2d>().swap(half);
		std::vector<Vector2d>().swap(next_half);

		entry_d_0 = Circuit::FindNearestPointPar(entry_spiral[entry_spiral.size() - 1], offset_0);
		exit_d_0 = Circuit::FindNearestPointPar(exit_spiral[exit_spiral.size() - 1], offset_0);

		Circuit::SelectOnePartOffset(offset_0, entry_d_0, exit_d_0, half);

		for (int i = 1; i < half.size(); i++)
		{
			entry_spiral.push_back(half[i]);
		}

	}

	void ToolpathGenerator::FermatsSpiralSmooth(std::vector<Vector2d> &contour, Vector2d input_entry_point, Vector2d input_exit_point)
	{
		ComputeOffsets(contour);

		entry_d_0 = Circuit::FindNearestPointPar(input_entry_point, offsets[0]);
		exit_d_0 = Circuit::FindNearestPointPar(input_exit_point, offsets[0]);

		for (int i = 0; i < offsets.size()-1; i++)
		{
			std::vector<Vector2d> current_half;
			std::vector<Vector2d> next_half;

			Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_0, current_half);
			for (int j = 0; j < current_half.size(); j++)
			{
				double d = Circuit::FindNearestPointPar(current_half[j], offsets[i + 1]);
				Vector2d v = Circuit::GetOnePointFromOffset(d, offsets[i + 1]);
				next_half.push_back(v);
			}

			Strip::StripDeformation(current_half, next_half);
			for (int j = 0; j < current_half.size(); j++)
			{
				if (i % 2 == 0)
					entry_spiral.push_back(current_half[j]);
				else
					exit_spiral.push_back(current_half[j]);
			}
			std::vector<Vector2d>().swap(current_half);
			std::vector<Vector2d>().swap(next_half);

			Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_0, current_half);

			for (int j = 0; j < current_half.size(); j++)
			{
				double d = Circuit::FindNearestPointPar(current_half[j], offsets[i + 1]);
				Vector2d v = Circuit::GetOnePointFromOffset(d, offsets[i + 1]);
				next_half.push_back(v);
			}

			Strip::StripDeformation(current_half, next_half);
			for (int j = 0; j < current_half.size(); j++)
			{
				if (i % 2 == 0)
					exit_spiral.push_back(current_half[j]);
				else
					entry_spiral.push_back(current_half[j]);
			}
			std::vector<Vector2d>().swap(current_half);
			std::vector<Vector2d>().swap(next_half);

			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offsets[i]);
			Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, offsets[i]);

			entry_d_0 = Circuit::FindNearestPointPar(entry_p_0, offsets[i + 1]);
			exit_d_0 = Circuit::FindNearestPointPar(exit_p_0, offsets[i + 1]);

		}

		std::vector<Vector2d> half_0;
		std::vector<Vector2d> half_1;

		Circuit::SelectOnePartOffset(offsets[offsets.size() - 1], entry_d_0, exit_d_0, half_0);
		Circuit::SelectOnePartOffset(offsets[offsets.size() - 1], exit_d_0, entry_d_0, half_1);

		if (Strip::GetTotalLength(half_0) > Strip::GetTotalLength(half_1))
		{

		}

		std::vector<Vector2d>().swap(half_0);
		std::vector<Vector2d>().swap(half_1);

	}

	void ToolpathGenerator::FermatSpiral(std::vector<Vector2d> &contour, Vector2d input_entry_point, Vector2d input_exit_point)
	{
		ComputeOffsets(contour);

		entry_d_0 = Circuit::FindNearestPointPar(input_entry_point, offsets[0]);
		exit_d_0 = Circuit::FindNearestPointPar(input_exit_point, offsets[0]);

		bool b0 = false;

		std::vector<Vector2d> offset;
		Circuit::GenerateOffset(offsets[offsets.size() - 1], toolpath_size / 2, offset);
		if (offset.size() == 0)
			b0=true;
		std::vector<Vector2d>().swap(offset);

		for (int i = 0; i < offsets.size()-1; i++)
		{

			if (i == offsets.size() - 2)
			{
				int dsad = 10;
			}

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

			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, offsets[i]);
			Vector2d entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, offsets[i]);
			Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, offsets[i]);
			Vector2d exit_p_1 = Circuit::GetOnePointFromOffset(exit_d_1, offsets[i]);

			if (i < offsets.size() - 1)
			{
				double entry_d_1_0 = Circuit::FindNearestPointPar(entry_p_1, offsets[i + 1]);
				double exit_d_1_0 = Circuit::FindNearestPointPar(exit_p_1, offsets[i + 1]);

				if (abs(entry_d_1_0 - exit_d_1_0) < 0.000001)
				{
					double t = DeltaDEuclideanDistance(exit_d_1_0, toolpath_size, i + 1);

					Vector2d v = Circuit::GetOnePointFromOffset(t, offsets[i + 1]);
					t = Circuit::FindNearestPointPar(v, offsets[i]);

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
				exit_p_1=Circuit::GetOnePointFromOffset(exit_d_1,offsets[i]);
				entry_p_1 = Circuit::GetOnePointFromOffset(entry_d_1, offsets[i]);
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
				double t = Circuit::Distance(v, offsets[i]);
				double t0 = Circuit::FindNearestPointPar(v, offsets[i]);
				
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
				turning_points_entry.push_back(entry_p_0);
				turning_points_exit.push_back(entry_p_1);
				turning_points_exit.push_back(exit_p_0);
				turning_points_entry.push_back(exit_p_1);

				if (i == offsets.size() - 1)
				{
					if (b0)
					{
						entry_spiral.push_back(entry_p_0);
						exit_spiral.push_back(exit_p_0);
						break;
						if (b1)
						{
							Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, exit_spiral);
							entry_spiral.push_back(entry_p_0);
						}
						else
						{
							Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, entry_spiral);
							exit_spiral.push_back(exit_p_0);
						}
					}
					else
					{
						if (b2)
						{
							Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, entry_spiral);
							Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, exit_spiral);
						}
						else
						{
							if (b1)
							{
								Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, exit_spiral);
								entry_spiral.push_back(entry_p_0);
							}
							else
							{
								Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, entry_spiral);
								exit_spiral.push_back(exit_p_0);
							}
						}
					}
				}
				else
				{
					Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, entry_spiral);
					Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, exit_spiral);
				}

				if (i + 1 < offsets.size())
				{
					entry_d_0 = Circuit::FindNearestPointPar(exit_spiral[exit_spiral.size() - 1], offsets[i + 1]);
					exit_d_0 = Circuit::FindNearestPointPar(entry_spiral[entry_spiral.size() - 1], offsets[i + 1]);
				}
			}
			else
			{
				turning_points_exit.push_back(entry_p_0);
				turning_points_entry.push_back(entry_p_1);
				turning_points_entry.push_back(exit_p_0);
				turning_points_exit.push_back(exit_p_1);

				if (i == offsets.size() - 1)
				{
					if (b0)
					{
						exit_spiral.push_back(entry_p_0);
						entry_spiral.push_back(exit_p_0);
						break;
						if (b1)
						{
							Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, entry_spiral);
							exit_spiral.push_back(entry_p_0);
						}
						else
						{
							Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, exit_spiral);
							entry_spiral.push_back(exit_p_0);
						}
					}
					else
					{
						if (b2)
						{
							Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, exit_spiral);
							Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, entry_spiral);
						}
						else
						{
							if (b1)
							{
								Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, entry_spiral);
								exit_spiral.push_back(entry_p_0);
							}
							else
							{
								Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, exit_spiral);
								entry_spiral.push_back(exit_p_0);
							}
						}
					}

				}
				else
				{
					Circuit::SelectOnePartOffset(offsets[i], entry_d_0, exit_d_1, exit_spiral);
					Circuit::SelectOnePartOffset(offsets[i], exit_d_0, entry_d_1, entry_spiral);
				}

				if (i + 1 < offsets.size())
				{
					entry_d_0 = Circuit::FindNearestPointPar(entry_spiral[entry_spiral.size() - 1], offsets[i + 1]);
					exit_d_0 = Circuit::FindNearestPointPar(exit_spiral[exit_spiral.size() - 1], offsets[i + 1]);
				}
			}

			if (!goon)
				break;
		}
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
			assert(false);
		}

		double returnd=DeltaDEuclideanDistance(d, distance, contour);
		std::vector<Vector2d>().swap(contour);
		return returnd;
	}

	double ToolpathGenerator::DeltaDEuclideanDistance(double d, double distance, std::vector<Vector2d> &contour)
	{
		Vector2d v=Circuit::GetOnePointFromOffset(d,contour);

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
				double delta_d = Circuit::FindNearestPointPar(vecs[i], contour);

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
				double delta_d = Circuit::FindNearestPointPar(vecs[i], contour);

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


}