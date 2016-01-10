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

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <iostream>


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


struct FaceInfo2
{
	FaceInfo2(){}
	int nesting_level;

	bool in_domain(){
		return nesting_level % 2 == 1;
	}
};

typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point                                                Point;



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

	struct TrunkNode
	{
		TrunkNode(int offset_id)
		{
			related_offset_id = offset_id;
		}
		int related_offset_id;
		//related_offset_id is not the trunknode id
		//trunk id is its index in trunknodes

		std::vector<double> connecting_leaf_nodes_points;
		std::vector<int> connecting_leaf_nodes_spiral_id;

		std::vector<double> connecting_trunk_nodes_points;
		std::vector<int> connecting_trunk_nodes_id;

		std::vector<double> connecting_points;
		std::vector<int> connecting_points_related_pathes_id;
		std::vector<int> connecting_points_pairs;

		//There are two kinds of connecting points, connecting leaf node and trunk node, connecting trunk node and trunk node.
		//Connecting_points includes all these two connecting points.
		//Each connecting point refers to a one_path in pathes. 

		void AllConnectingPoints()
		{
			for (int i = 0; i < connecting_leaf_nodes_points.size(); i++)
			{
				connecting_points.push_back(connecting_leaf_nodes_points[i]);
			}
			for (int i = 0; i < connecting_trunk_nodes_points.size(); i++)
			{
				//connecting_points.push_back(connecting_trunk_nodes_points[i]);
			}

			for (int i = 0; i < connecting_points.size(); i++)
			{
				connecting_points_related_pathes_id.push_back(-1);
			}
		}
	};


	class ToolpathGenerator
	{
	protected:
		TMeshProcessor* m_render;
		
		Polygon_with_holes contours;
		
		double contour_size;
		
		std::vector<std::vector<Vector2d>> offsets;
		std::vector<std::vector<std::vector<Vector2d>>> offsetses;

		std::vector<std::vector<std::vector<Vector2d>>> offsets_parts;

		std::vector<int> offset_graph;

		std::vector<int> offset_degree;
		std::vector<std::vector<int>> decompose_offset;

		bool use_save_offset_file;
		std::string offset_file;

		std::vector<Vector2d> one_single_path;

		ImageSpace image_space;

		double toolpath_size;
		double input_toolpath_size;
		int input_int_1;
		int input_int_2;
		int line_width;
		int point_size;

		int work_model = 0;

		std::string load_path;

		int debug_int_0;
		int debug_int_1;
		int debug_int_2;

		std::vector<TrunkNode> trunk_nodes;

	public:
		
		bool draw_pixels;
		bool draw_aixs;
		bool draw_contour;
		bool draw_offsets;
		bool draw_turning_points;
		bool draw_entry_exit_spiral;
		bool draw_spiral;

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

		std::vector<Vector2d> debug_points;
		std::vector<std::vector<Vector2d>> debug_lines;

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

		void FermatsSpiralTrick(std::vector<std::vector<Vector2d>> &local_offsets, 
			Vector2d input_entry_point, Vector2d input_exit_point, Vector2d &output_entry_point, Vector2d &output_exit_point);

		void GenerateZigzag();
		void GenerateZigzagForArbitraryShape();
		void ArchinedeanSpiral(std::vector<Vector2d> &contour);
		void ArchinedeanSpiralSmooth(std::vector<Vector2d> &contour);
		void ArchinedeanSpiralTrick(std::vector<Vector2d> &contour);

		double ComputeNextTurningPoint(double d, double distance, int offset_index);

		void ContourSmoothing(std::vector<Vector2d> &contour);
		void DirectlyContourSmoothing(std::vector<Vector2d> &contour);

		void PolygonSmoothing();
		void OptimalDirection();
		void DirectlyPolygonSmoothing();


		//filling algorithm

		void BuildOffsetGraph();
		void FillingAlgorithmBasedOnOffsets();

		void Input_Offsetses(std::string path);
		void TurnPath();


		//output
		void Output_Obj(std::string path);
		void Output_Obj_1(std::string path);
		void OutputPath(std::vector<Vector2d> &vecs, std::string path);
		void OutputPathTwoCircles();

		void Output_tree(std::string path);
		void Output_tree(std::vector<int> &nodes, std::vector<int> &edges, std::string path);
		void Output_Offsetses(std::string path);

		void OutputPath_Hollow(std::string path);



		//triangle
		void mark_domains(CDT& ct, CDT::Face_handle start,int index, std::list<CDT::Edge>& border);
		void mark_domains(CDT& cdt);
		void insert_polygon(CDT& cdt, const Polygon_2& polygon, std::vector<int> &indexInt);

	};
} 