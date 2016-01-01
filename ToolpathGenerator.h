#pragma once

#include "MeshObject.h"
#include "MathHelper.h"
#include "OpenMeshHelper.h"
#include "MyMesh.h"

#include <iostream>
#include <fstream>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>

#include<boost/shared_ptr.hpp>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_with_holes_2.h>
#include<CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/barycenter.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/intersections.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>

//#include "Circuit.h"

typedef CGAL::Simple_cartesian<double> KC;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_without_intersections_2<KC> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>  SDG2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                    Point_2;
typedef K::Line_2                     Line_2;
typedef K::Segment_2                  Segment_2;
typedef CGAL::Polygon_2<K>            Polygon_2;

typedef CGAL::Straight_skeleton_2<K>  Ss;
typedef Ss::Vertex_const_handle     Vertex_const_handle;
typedef Ss::Halfedge_const_handle   Halfedge_const_handle;
typedef Ss::Halfedge_const_iterator Halfedge_const_iterator;

typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes;
typedef CGAL::Straight_skeleton_2<K>  Ss;
typedef boost::shared_ptr<Polygon_2> PolygonPtr;
typedef boost::shared_ptr<Ss> SsPtr;
typedef std::vector<PolygonPtr> PolygonPtrVector;

#define MAXDOUBLE 1000000000.0;

namespace hpcg {

	class TMeshProcessor;

	struct Pixel
	{
		Vector2d center;
		bool inside;
		bool filled;

		void reset()
		{
			filled = false;
			inside = false;
		}

		Pixel(Vector2d d)
		{
			center = d;
			reset();
		}
	};

	struct ImageSpace
	{
		std::vector<std::vector<Pixel>> pixels;
		double pixel_size;
		double toolpath_size;

		bool check(int i,int j)
		{
			return (i >= 0 && i < pixels.size() && j >= 0 && j < pixels[0].size());
		}

		void Reset()
		{
			for (int i = 0; i < pixels.size(); i++)
			{
				for (int j = 0; j < pixels[i].size(); j++)
				{
					pixels[i][j].reset();
				}
			}
		}
	};

	class ToolpathGenerator
	{
	protected:
		TMeshProcessor* m_render;
		
		Polygon_with_holes contours;
		std::vector<std::vector<Vector2d>> offsets;
		std::vector<std::vector<std::vector<Vector2d>>> offsetses;
		std::vector<int> offset_graph;
		std::vector<int> offset_degree;
		std::vector<std::vector<int>> decompose_offset;

		ImageSpace image_space;

		double toolpath_size;
		int input_int_1;
		int input_int_2;
		int line_width;
		int point_size;

		int work_model = 0;

		std::string load_path;

		int debug_int_0;
		int debug_int_1;

	public:
		
		bool draw_pixels;
		bool draw_aixs;
		bool draw_contour;
		bool draw_offsets;
		bool draw_spiral;
		bool draw_turning_points;
		bool draw_voronoi;
		bool draw_medial_axis;
		bool draw_minimal_points;
		bool draw_entry_exit_spiral;
		bool draw_entry_exit_points;
		bool draw_polygons_entry_exit;
		bool draw_cutting_points;
		bool draw_inner_concave_points;
		bool draw_maximal_points;
		bool draw_save_critical_points;

		int smooth_number = 0;

		SDG2 sdg;

		std::vector<Vector2d> entry_spiral;
		std::vector<Vector2d> exit_spiral;

		std::vector<std::vector<Vector2d>> entry_spirals;
		std::vector<std::vector<Vector2d>> exit_spirals;

		std::vector<std::vector<Vector2d>> pathes;
		std::vector<std::vector<Vector2d>> pathes_temp;

		std::vector<Vector2d> turning_points_entry;
		std::vector<Vector2d> turning_points_exit;

		std::vector<Vector2d> turning_points_entry_temp;
		std::vector<Vector2d> turning_points_exit_temp;

		double entry_d_0;
		double exit_d_0;

		std::vector<Vector2d> aaaa;



		ToolpathGenerator();
		void Init(TMeshProcessor* render);
		void LoadContour();

		void BuildeImageSpace();
		void StepDebug();

		//rending
		void DrawPixel(int i,int j);
		void Rendering();

		//offset fermat spiral
		void FermatSpiral(std::vector<Vector2d> &contour, Vector2d input_entry_point, Vector2d input_exit_point);

		void FermatsSpiralTrick(std::vector<std::vector<Vector2d>> &local_offsets, Vector2d input_entry_point, Vector2d input_exit_point);

		void GenerateZigzag();
		void GenerateZigzagForCircle();
		void ArchinedeanSpiral(std::vector<Vector2d> &contour);
		void ArchinedeanSpiralSmooth(std::vector<Vector2d> &contour);
		void ArchinedeanSpiralTrick(std::vector<Vector2d> &contour);

		double ComputeNextTurningPoint(double d, double distance, int offset_index);

		void PolygonSmoothing();
		void OutputPath(std::vector<Vector2d> &vecs, std::string path);
		void OutputPath(std::string path);
		void OutputPathTwoCircles();

		//filling algorithm
		void FillingAlgorithmBasedOnOffsets();

		void Output_tree(std::string path);
		void Output_tree(std::vector<int> &nodes, std::vector<int> &edges, std::string path);

	};
} 