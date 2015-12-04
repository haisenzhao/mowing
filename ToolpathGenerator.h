#pragma once

#include "MeshObject.h"
#include "MathHelper.h"
#include "OpenMeshHelper.h"
#include "MyMesh.h"

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
		bool boundary;
		bool centeraxis;

		int filled_spiral_index;
		bool filled;
		bool filled_entry_exit;//true:entry; false:exit;
		double distance;
		double distance_entry;
		double distance_exit;
		
		double save_distance;
		int save_related_boundary_index;
		int save_related_boundary_x;
		int save_related_boundary_y;

		bool bool_t;
		bool bool_t_entry_exit;
		bool bool_t_entry;
		bool bool_t_exit;

		int related_boundary_index;
		int related_boundary_index_entry;
		int related_boundary_index_exit;

		int related_boundary_x;
		int related_boundary_y;

		int related_boundary_x_entry;
		int related_boundary_y_entry;

		int related_boundary_x_exit;
		int related_boundary_y_exit;

		int boundary_index;


		int failure_nb;

		void reset()
		{
			boundary = false;
			bool_t = false;
			centeraxis = false;
			boundary_index = -1;
			distance = MAXDOUBLE;
			distance_entry = MAXDOUBLE;
			distance_exit = MAXDOUBLE;

			related_boundary_index_entry = -1;
			related_boundary_index_exit = -1;

			related_boundary_index = -1;
			related_boundary_x = -1;
			related_boundary_y = -1;
			filled_spiral_index = -1;

			bool_t_entry = false;
			bool_t_exit = false;

			failure_nb = 0;
		}

		Pixel(Vector2d d)
		{
			center = d;
			filled = false;
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

		void CheckBoundaryForOnePixel(int i, int j)
		{
			if (pixels[i][j].inside)
			{
				if (pixels[i][j].filled)
				{
					for (int a = -1; a <= 1; a++)
					{
						for (int b = -1; b <= 1; b++)
						{
							if (check(a + i, b + j) && !(a == 0 && b == 0) && a*b == 0)
							{
								if (!pixels[a + i][b + j].filled)
								{
									pixels[i][j].boundary = true;
								}
							}
						}
					}
				}
			}
			else
			{
				for (int a = -1; a <= 1; a++)
				{
					for (int b = -1; b <= 1; b++)
					{
						if (check(a + i, b + j) && !(a == 0 && b == 0) && a*b == 0)
						{
							if (pixels[a + i][b + j].inside)
							{
								pixels[i][j].boundary = true;
							}
						}
					}
				}
			}
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

		void CheckBoundary()
		{
			for (int i = 0; i < pixels.size(); i++)
			{
				for (int j = 0; j < pixels[i].size(); j++)
				{
					CheckBoundaryForOnePixel(i, j);
				}
			}
		}

		void UpdateDistanceFromBoundary()
		{
			for (int i = 0; i < pixels.size(); i++)
			{
				for (int j = 0; j < pixels[i].size(); j++)
				{
					if (pixels[i][j].boundary)
					{
						for (int k = i - toolpath_size / pixel_size; k < i + toolpath_size / pixel_size; k++)
						{
							for (int l = j - toolpath_size / pixel_size; l < j + toolpath_size / pixel_size; l++)
							{
								if (k == 27 && l == 5)
								{
									int dsd = 10;
								}

								if (check(i, j) && check(k, l))
								{
									if (pixels[k][l].distance>PixelsDistance(i, j, k, l))
									{
										pixels[k][l].distance = PixelsDistance(i, j, k, l);
										pixels[k][l].distance_entry = pixels[k][l].distance;
										pixels[k][l].distance_exit = pixels[k][l].distance;
										pixels[k][l].related_boundary_x = i;
										pixels[k][l].related_boundary_y = j;
										pixels[k][l].related_boundary_x_entry = i;
										pixels[k][l].related_boundary_y_entry = j;
										pixels[k][l].related_boundary_x_exit = i;
										pixels[k][l].related_boundary_y_exit = j;
									}
								}
							}
						}
					}
				}
			}
		}

		bool PixelsNeighbour(int i0, int j0, int i1, int j1)
		{
			return abs(i0 - i1) <= 1 && abs(j0 - j1) <= 1;
		}

		double PixelsDistance(int i0, int j0, int i1, int j1)
		{
			if (check(i0, j0) && check(i1, j1))
			{
				return sqrt((double)CGAL::squared_distance(Point_2(pixels[i0][j0].center[0], pixels[i0][j0].center[1]),
					Point_2(pixels[i1][j1].center[0], pixels[i1][j1].center[1])));
			}
			else
			{
				assert(false);
				return -1.0;
			}
		}

		void FindClosestPoint(double x, double y, int &index_x, int &index_y)
		{
			double min_d = MAXDOUBLE;

			for (int i = (int)((x - pixels[0][0].center[0]) / pixel_size) - 3; i <= (int)((x - pixels[0][0].center[0]) / pixel_size) + 3; i++)
			{
				for (int j = (int)((y - pixels[0][0].center[1]) / pixel_size) - 3; j <= (int)((y - pixels[0][0].center[1]) / pixel_size) + 3; j++)
				{
					if (check(i, j))
					{
						double l = sqrt((double)CGAL::squared_distance(Point_2(x, y), Point_2(pixels[i][j].center[0], pixels[i][j].center[1])));
						if (l < min_d)
						{
							index_x = i;
							index_y = j;
							min_d = l;
						}
					}
				}
			}
		}

	};

	struct PixelIndex
	{
		int x;
		int y;
		PixelIndex(int a,int b)
		{
			x = a;
			y = b;
		}
	};

	class ToolpathGenerator
	{
	protected:
		TMeshProcessor* m_render;
		
		Polygon_with_holes contours;
		std::vector<std::vector<Vector2d>> toolpath;
		ImageSpace image_space;
		std::vector<PixelIndex> spiral;
		std::vector<PixelIndex> unfilled_region_boundary;
		std::vector<std::vector<PixelIndex>> spirals;

		double toolpath_size;
		int spiral_affect_index;
		int input_int_1;
		int input_int_2;
		bool spiral_direction;

		std::vector<PixelIndex> spiral_entry;
		std::vector<PixelIndex> spiral_exit;

		int linewidth;

	public:
		
		bool draw_boundary;
		bool draw_pixels;
		bool draw_aixs;
		bool draw_contour;
		bool draw_offset;
		bool draw_spiral;

		SsPtr iss;

		ToolpathGenerator();
		void Init(TMeshProcessor* render);
		void LoadContour();
		void GenerateToolpath();
		void GenerateToolpathBasedonOffset();
		double SpiralPathDistance(int i, int j);

		void BuildeImageSpace();
		bool FindStartPixel(int &index_x,int &index_y);
		void InsertPixeltoSpiral(PixelIndex pi);
		
		void SearchUnfilledRegionBoundary(int start_x,int start_y);
		void SearchOnePath(bool b);
		void GenerateToolpathBasedonImageSpace();

		void StepDebug();

		void DrawPixel(int i,int j);
		void Rendering();

		//Fermat spiral
		int entry_point_x, entry_point_y;
		int exit_point_x, exit_point_y;

		int spiral_entry_affect_index;
		int spiral_exit_affect_index;

		bool first_two_step_bool;

		void FermatSpiral();
		void ChooseEntryExitPoints();
		void GetOnePointFromOffset(double d, int &index_x,int &index_y);
		
		void InsertPixeltoSpiralEntry(PixelIndex pi);
		void InsertPixeltoSpiralExit(PixelIndex pi);
		void SearchOneStepforSpiralEntry(bool b);
		void SearchOneStepforSpiralExit(bool b);

		void ConnectTwoPixels(PixelIndex p1, PixelIndex p2, std::vector<PixelIndex> &insertpixels);

		//offset fermat spiral
		void OffsetBasedFermatSpiral();
		void OffsetBasedFermatSpiral1();

		Vector2d entry_point;
		Vector2d exit_point;

		std::vector<Vector2d> entry_points;
		std::vector<Vector2d> exit_points;
		
		std::vector<Vector2d> entry_spiral;
		std::vector<Vector2d> exit_spiral;

		std::vector<Vector2d> aaaa;
		std::vector<Vector2d> bbbb;
		std::vector<Vector2d> cccc;
		std::vector<Vector2d> dddd;

		std::vector<Vector2d> contour_aaa;

		double entry_d_0;
		double exit_d_0;

		double OneOffsetLength(int offset_index);
		double FindNearestPoint(int offset_index, Vector2d v);
		double FindNearestPoint(Vector2d v, std::vector<Vector2d> &contour);

		double DeltaDGeodesicDistance(double d, double distance, std::vector<Vector2d> &contour);
		double DeltaDEuclideanDistance(double d, double distance, std::vector<Vector2d> &contour);
		double DeltaDEuclideanDistance(double d, double distance, int offset_index);
		double ComputeNextTurningPoint(double d, double distance, int offset_index);

		void GetOnePointFromOffset(int offset_index, double d, Vector2d &v);
		void GetOnePointFromOffset(std::vector<Vector2d> &contour, double d, Vector2d &v);
		
		void SelectOnePartOffset(int offset_index, double d0, double d1, std::vector<Vector2d> &vecs);
		void SelectOnePartOffset(std::vector<Vector2d> &contour, double d0, double d1, std::vector<Vector2d> &vecs);
		void GenerateOffset(bool direction, std::vector<Vector2d> &contour, double lOffset, std::vector<Vector2d> &offset);

		void GenerateOffset(bool direction, int offset_index, double lOffset, std::vector<Vector2d> &offset);

		bool OneStep(std::vector<Vector2d> &outside_contour, double &start_d, double &end_d, std::vector<Vector2d> &entry_spiral_points, std::vector<Vector2d> &exit_spiral_points);
		bool OneStep1(std::vector<Vector2d> &outside_contour, double &start_d, double &end_d, std::vector<Vector2d> &entry_spiral_points, std::vector<Vector2d> &exit_spiral_points);

		double ComputeDistancePointContour(std::vector<Vector2d> &contour, Vector2d v);

		bool  TwoContourIntersection(std::vector<Vector2d> &contour0, std::vector<Vector2d> &contour1);

		double MinimalDistance(std::vector<Vector2d> &vecs, Vector2d &v);

		void CreateMAT();

		void Zigzag();
		void ArchinedeanSpiral();

		void PolygonSmoothing();
	};
} 