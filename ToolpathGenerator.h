#pragma once

#include "MeshObject.h"
#include "MathHelper.h"
#include "OpenMeshHelper.h"
#include "MyMesh.h"

#include <iostream>
#include <fstream>
#include <cassert>

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

	struct MedialAxis
	{
		//std::vector<Vector2d> 
	};

	class ToolpathGenerator
	{
	protected:
		TMeshProcessor* m_render;
		
		Polygon_with_holes contours;
		std::vector<std::vector<Vector2d>> offsets;
		ImageSpace image_space;

		double toolpath_size;
		int input_int_1;
		int input_int_2;
		int linewidth;

	public:
		
		bool draw_boundary;
		bool draw_pixels;
		bool draw_aixs;
		bool draw_contour;
		bool draw_offset;
		bool draw_spiral;

		SDG2 sdg;

		std::vector<Vector2d> entry_spiral;
		std::vector<Vector2d> exit_spiral;
		std::vector<Vector2d> cccc;
		std::vector<Vector2d> dddd;
		double entry_d_0;
		double exit_d_0;

		std::vector<Vector2d> mats;

		ToolpathGenerator();
		void Init(TMeshProcessor* render);
		void LoadContour();
		void ComputeOffsets();

		void BuildeImageSpace();
		void StepDebug();

		//rending
		void DrawPixel(int i,int j);
		void Rendering();

		//offset fermat spiral
		void GenerateFermatSpiral();
		void GenerateZigzag();
		void ArchinedeanSpiral();

		double FindNearestPointPar(Vector2d v, int offset_index);
		double FindNearestPointPar(Vector2d v, std::vector<Vector2d> &contour);

		double DeltaDGeodesicDistance(double d, double distance, std::vector<Vector2d> &contour);
		double DeltaDEuclideanDistance(double d, double distance, std::vector<Vector2d> &contour);
		double DeltaDEuclideanDistance(double d, double distance, int offset_index);
		double ComputeNextTurningPoint(double d, double distance, int offset_index);

		Vector2d GetOnePointFromOffset(double d, int offset_index);
		Vector2d GetOnePointFromOffset(double d,std::vector<Vector2d> &contour);
		
		void SelectOnePartOffset(int offset_index, double d0, double d1, std::vector<Vector2d> &vecs);
		void SelectOnePartOffset(std::vector<Vector2d> &contour, double d0, double d1, std::vector<Vector2d> &vecs);
		
		void GenerateOffset(bool direction, std::vector<Vector2d> &contour, double lOffset, std::vector<Vector2d> &offset);
		void GenerateOffset(bool direction, int offset_index, double lOffset, std::vector<Vector2d> &offset);

		double MinimalDistance(std::vector<Vector2d> &vecs, Vector2d &v);
		double MinimalDistanceSegments(std::vector<Vector2d> &segments, Vector2d &v);

		double MinimalDistanceContours(Vector2d &v);

		void CreateDelaunayGraphs();

		void PolygonSmoothing();

		void OutputPath(std::vector<Vector2d> &vecs, std::string path);

		bool CheckInsideContours(Vector2d &v);

		bool CheckOnContours(Vector2d &v);

		void RelatedPointsOnContours(Vector2d &v, double distace, std::vector<Vector2d> &vecs);

	};
} 