#include "t8dg.h"
#include "t8dg_quadrature.h"
#include "t8dg_vertexset.h"

/** Additionally to the vertices save the quadrature weights*/
struct t8dg_quadrature
{
  int                 dim;
  int                 number_of_quadrature_points;     /**< Number of element quadrature points*/
  t8dg_vertexset_t   *vertices;                       /**< LGL quadrature vertices*/
  sc_array_t         *weights;                          /**< LGL quadrature weights*/
};

t8dg_quadrature_t  *
t8dg_quadrature_new (t8dg_vertexset_t * vertexset)
{
  T8DG_CHECK_ABORT (t8dg_vertexset_get_dim (vertexset) == 1 && t8dg_vertexset_get_type (vertexset) == T8DG_LGL, "Not yet implemented");

  t8dg_quadrature_t  *quadrature = T8DG_ALLOC (t8dg_quadrature_t, 1);

  quadrature->dim = t8dg_vertexset_get_dim (vertexset);
  quadrature->vertices = vertexset;
  t8dg_vertexset_ref (vertexset);
  quadrature->weights = sc_array_new_count (sizeof (double), t8dg_vertexset_get_num_element_vertices (vertexset));

  double             *weights;
  weights = (double *) t8_sc_array_index_locidx (quadrature->weights, 0);

  /* sum of element weights = 1 = Vol([0,1]) */
  switch (t8dg_vertexset_get_num_element_vertices (vertexset)) {
  case (1):
    weights[0] = 1;
    break;
  case (2):
    weights[0] = 0.5;
    weights[1] = 0.5;
    break;
  case (3):
    weights[0] = 1. / 6;
    weights[1] = 4. / 6;
    weights[2] = 1. / 6;
    break;
  case (4):
    weights[0] = 1. / 12;
    weights[1] = 5. / 12;
    weights[2] = 5. / 12;
    weights[3] = 1. / 12;
    break;
  default:
    T8DG_ABORT ("Not yet implemented!\n");
  }
  return quadrature;
}

void
t8dg_quadrature_destroy (t8dg_quadrature_t ** pquadrature)
{
  t8dg_quadrature_t  *quadrature = *pquadrature;
  t8dg_vertexset_unref (quadrature->vertices);
  quadrature->vertices = NULL;
  sc_array_destroy (quadrature->weights);
  quadrature->weights = NULL;
  quadrature->dim = -1;
  quadrature->number_of_quadrature_points = -1;
  T8_FREE (quadrature);
  *pquadrature = NULL;
}

int
t8dg_quadrature_get_num_faces (const t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);
  return t8dg_vertexset_get_num_faces (quadrature->vertices);
}

t8dg_quad_idx_t
t8dg_quadrature_get_num_element_vertices (const t8dg_quadrature_t * quadrature)
{
  T8DG_ASSERT (quadrature != NULL);
  return t8dg_vertexset_get_num_element_vertices (quadrature->vertices);
}

t8dg_quad_idx_t
t8dg_quadrature_get_num_face_vertices (const t8dg_quadrature_t * quadrature, const int iface)
{
  T8DG_ASSERT (quadrature != NULL);
  return t8dg_vertexset_get_num_face_vertices (quadrature->vertices, iface);
}

void
t8dg_quadrature_get_element_vertex (double vertex[3], const t8dg_quadrature_t * quadrature, const t8dg_quad_idx_t iquad)
{
  T8DG_ASSERT (quadrature != NULL);
  t8dg_vertexset_get_vertex (vertex, quadrature->vertices, iquad);
}

double
t8dg_quadrature_get_element_weight (const t8dg_quadrature_t * quadrature, const t8dg_quad_idx_t iquad)
{
  T8DG_ASSERT (quadrature != NULL);
  return *(double *) t8dg_sc_array_index_quadidx (quadrature->weights, iquad);
}

void
t8dg_quadrature_get_face_vertex (double vertex[3], const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t ifacequad)
{
  T8DG_CHECK_ABORT (t8dg_vertexset_get_type (quadrature->vertices) == T8DG_LGL, "Not yet implemented");
  T8DG_ASSERT (quadrature != NULL);
  t8dg_quad_idx_t     ielemquad;
  ielemquad = t8dg_vertexset_get_LGL_facevertex_element_index (quadrature->vertices, iface, ifacequad);
  t8dg_vertexset_get_vertex (vertex, quadrature->vertices, ielemquad);
}

double
t8dg_quadrature_get_face_weight (const t8dg_quadrature_t * quadrature, const int iface, const t8dg_quad_idx_t ifacequad)
{
  T8DG_CHECK_ABORT (quadrature->dim == 1, "Not yet implemented");
  return 1;
}

t8dg_vertexset_type_t
t8dg_quadrature_get_type (const t8dg_quadrature_t * quadrature)
{
  return t8dg_vertexset_get_type (quadrature->vertices);
}

t8dg_vertexset_t   *
t8dg_quadrature_get_vertexset (const t8dg_quadrature_t * quadrature)
{
  return quadrature->vertices;
}