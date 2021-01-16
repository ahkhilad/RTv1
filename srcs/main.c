/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/01/04 21:42:11 by ahkhilad          #+#    #+#             */
/*   Updated: 2021/01/16 19:41:52 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

int ft_is_zero(float value){
	return (value < 0.5f  && value > -0.001f);
}

void ft_print_vect(t_vec v){

	printf("vec : x : %.2f\t y : %.2f\t z : %.2f\n", v.x, v.y, v.z);
}

static int		ft_min_ray(float t1, float t2, float *t)
{
	if (((t1 < t2 || t2 < 0.001) && t1 > 0.1) && (t1 < *t))
	{
		*t = t1;
		return (1);
	}
	else if (((t2 < t1 || t1 < 0.001) && t2 > 0.1) && (t2 < *t))
	{
		*t = t2;
		return (1);
	}
	else
		return (0);
}

int sphere_intersect(t_object *sphere, t_ray *ray, float *tmin)
{
	t_vec x;
	float a, b, c, delta;
	float t1, t2;

	x = ft_vectorsub(ray->source, sphere->pos);
	a = ft_dotproduct(ray->direction, ray->direction);
	b = 2.0 * ft_dotproduct(ray->direction, x);
	c = ft_dotproduct(x, x) - sphere->radius * sphere->radius;
	delta = (b * b) - (4.0 * a * c);
	if (delta < 0)
		return (0);
	delta = sqrtf(delta);
	t1 = (-b + delta) / (2 * a);
	t2 = (-b - delta) / (2 * a);
	return (ft_min_ray(t1, t2, tmin));
}

int		plane_intersect(t_object *plane, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	a;
	float	b;
	float	t1;

	x = ft_vectorsub(ray->source, plane->pos);
	a = -1.0 * ft_dotproduct(x, plane->normal);
	b = ft_dotproduct(ray->direction, plane->normal);
	if (fabs(b) <= 1e-6)
		return (0);
	t1 = a / b;
	if (t1 > 1e-2 && t1 < *tmin)
	{
		*tmin = t1;
		return (1);
	}
	return (0);
}

static int	ft_cylindre_cap(t_ray *ray, float *t, t_vec pos, t_vec n)
{
	float	ddn;
	float	t1;
	t_vec	dist;

	ddn = ft_dotproduct(ray->direction, n);
	if (ddn <= 1.0e-6 && ddn >= -1.0e-6)
		return (0);
	dist = ft_vectorsub(pos, ray->source);
	t1 = (ft_dotproduct(dist, n)) / ddn;
	// if (t1 < *t && t1 > 0) // more stuff needs to be re_do
	if (t1 >= 0.0)
	{
		*t = t1;
		if (ddn >= 1e-6)
		 	return (2);
		return (1);
	}
	return (0);
}

static int	ft_cylinder_cap1(t_object *cylinder, t_ray *ray, float *t, float m1)
{
	if (m1 < -cylinder->height)
		return (0);
	if (ft_cylindre_cap(ray, t, ft_vectormulti(cylinder->axis, -cylinder->height), ft_vectormulti(cylinder->axis, -1.0)))
		return (-1);
	return (0);
}

int		cylinder_intersect(t_object *cylinder, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	a;
	float	b;
	float	c;
	float	delta;
	float	t1;
	float	t2;
	float	m0;
	float	m1;

	// printf("%.3f", cylinder->height);
	// exit(0);
	
	x = ft_vectorsub(ray->source, cylinder->pos);
	//printf("x : %.2f\t y : %.2f \t z = %.2f\n", ray->direction.x, ray->direction.y, ray->direction.z);
	//x = ft_vectorsub( cylinder->pos, ray->source);
	a = ft_dotproduct(ray->direction, ray->direction) - pow(ft_dotproduct(ray->direction, cylinder->axis), 2.0);
	b = 2.0 * (ft_dotproduct(ray->direction, x) - (ft_dotproduct(ray->direction, cylinder->axis) * ft_dotproduct(x, cylinder->axis)));
	c = ft_dotproduct(x, x) - pow(ft_dotproduct(x, cylinder->axis), 2.0) - (cylinder->radius * cylinder->radius);
	delta = (b * b) - (4.0 * a * c);
	//printf("a : %.2f\t b : %.2f \t c = %.2f\n", a, b, c);
	if (delta < 0)
		return (0);
	delta = sqrt(delta);
	t1 = (-b + delta) / (2 * a);
	t2 = (-b - delta) / (2 * a);
	if (cylinder->height <= 0)
		return (ft_min_ray(t1, t2, tmin));
	if (t1 > t2)
	{
		m0 = t1;
		t1 = t2;
		t2 = m0;
	}
	m0 = ft_dotproduct(ray->direction, cylinder->axis) * t1;
	m0 += ft_dotproduct(x, cylinder->axis);
	m1 = ft_dotproduct(ray->direction, cylinder->axis) * t2;
	m1 += ft_dotproduct(x, cylinder->axis);
    //printf("m0 : %.2f\t  m1 : %.2f\n", m0, m1);
	if (m0 < -cylinder->height)
		return (ft_cylinder_cap1(cylinder, ray, tmin, m1));
	else if (m0 >= -cylinder->height && m0 <= cylinder->height)
	{
		if (t1 < 1e-6 || t1 > *tmin)
			return (0);
		*tmin = t1;
		return (1);
	}
	else if (m0 > cylinder->height)
	{
		if (m1 > cylinder->height)
			return (0);
		if (ft_cylindre_cap(ray, tmin, ft_vectormulti(cylinder->axis, cylinder->height), cylinder->axis))
			return (2);
	}
	return (0);
}

float			ft_get_cone_limits(float t, t_ray *r, t_vec va, t_vec v)
{
	t_vec	q;
	float	m1;
	float	m2;

	if (t < 0)
		return (INFINITY);
	q = ft_vectoradd(r->source, ft_vectormulti(r->direction, t));
	m1 = ft_dotproduct(va, q);
	m2 = ft_dotproduct(va, ft_vectorsub(q, v));
	if (m1 < 0.0 && m2 > 0)
		return (t);
	return (INFINITY);
}

typedef struct	s_delta
{
	float	a;
	float	b;
	float	c;
	float	delta;
}				t_delta;

typedef struct	s_aux_cone
{
	t_delta	d;
	t_vec	pa;
	t_vec	va;
	t_vec	x;
	t_vec	a;
	t_vec	b;
	float	anglesin;
	float	t1;
	float	t2;
	float	m;
	float	anglecos;
}				t_aux_cone;

typedef struct	s_aux_parallelogram
{
	t_vec	e01;
	t_vec	e03;
	t_vec	e23;
	t_vec	e21;
	t_vec	p;
	t_vec	p2;
	t_vec	t;
	t_vec	t2;
	t_vec	q;
	t_vec	q2;
	t_vec	n;
	float	det;
	float	det2;
	float	invdet;
	float	invdet2;
	float	t1;
	float	a;
	float	b;
	float	a2;
	float	b2;
}				t_aux_parallelogram;

static void		ft_aux_cone_init(t_aux_cone *a, t_object *c, t_ray *r)
{
	a->anglecos = pow(cos(c->angle), 2.0);
	a->anglesin = pow(sin(c->angle), 2.0);
	a->pa = (t_vec){0.0, 0.0, 0.0};
	a->va = ft_normalize(ft_vectorsub(a->pa, c->axis));
	a->x = ft_vectorsub(r->source, a->pa);
	a->a = ft_vectorsub(r->direction, ft_vectormulti(a->va, ft_dotproduct(r->direction, a->va)));
	a->b = ft_vectorsub(a->x, ft_vectormulti(a->va, ft_dotproduct(a->x, a->va)));
	a->d.a = a->anglecos * ft_dotproduct(a->a, a->a) -
		a->anglesin * ft_dotproduct(r->direction, a->va) * ft_dotproduct(r->direction, a->va);
	a->d.b = 2.0 * a->anglecos * ft_dotproduct(a->a, a->b) -
		2.0 * a->anglesin * ft_dotproduct(r->direction, a->va) * ft_dotproduct(a->x,
				a->va);
	a->d.c = a->anglecos * ft_dotproduct(a->b, a->b) -
		a->anglesin * ft_dotproduct(a->x, a->va) * ft_dotproduct(a->x, a->va);
	a->d.delta = a->d.b * a->d.b - 4.0 * a->d.a * a->d.c;
}

int				cone_intersect(t_object *c, t_ray *r, float *t)
{
	t_aux_cone	a;

	ft_aux_cone_init(&a, c, r);
	if (a.d.delta < 0.00000001)
		return (0);
	a.d.delta = sqrt(a.d.delta);
	if (c->height <= 0)
	{
		return (ft_min_ray((-a.d.b + a.d.delta) / (2.0 * a.d.a),
			(-a.d.b - a.d.delta) / (2.0 * a.d.a), t));
	}
	a.t1 = ft_get_cone_limits((-a.d.b + a.d.delta) / (2.0 * a.d.a),
			r, a.va, c->axis);
	a.t2 = ft_get_cone_limits((-a.d.b - a.d.delta) / (2.0 * a.d.a),
			r, a.va, c->axis);
	if (!ft_min_ray(a.t1, a.t2, t))
		return (0);
	a.m = ft_dotproduct(r->direction, c->axis) * (*t);
	a.m += ft_dotproduct(r->source, c->axis);
	if (a.m > c->height / 2.0 - 5.0e-1)
		return (2);
	return (1);
}

/*
int		cone_intersect(t_object *cone, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	k;
	float	a;
	float	b;
	float	c;
	float	delta;
	float	m0;
	float	m1;
	float	t1;
	float	t2;

	x = ft_vectorsub(ray->source, cone->pos);
	//x = ft_vectorsub( cone->pos, ray->source);
	k = tanf(deg_to_rad(cone->angle) / 2.0);
	a = ft_dotproduct(ray->direction, ray->direction) - (1.0 + (k * k)) * powf(ft_dotproduct(ray->direction, cone->axis), 2.0);
	b = 2.0 * (ft_dotproduct(ray->direction, x) - ((1.0 + (k * k)) * ft_dotproduct(ray->direction, cone->axis) * ft_dotproduct(x, cone->axis)));
	c = ft_dotproduct(x, x) - (1.0 + (k * k)) * powf(ft_dotproduct(x, cone->axis), 2.0);
	delta = (b * b) - (4.0 * a * c);
	if (delta < 0)
		return (0);
	delta = sqrtf(delta);
	t1 = (-b + delta) / (2 * a);
	t2 = (-b - delta) / (2 * a);
	return (ft_min_ray(t1, t2, tmin));
	
}*/

int		ellipsoid_intersect(t_object *ellipsoid, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	r;
	float	a1;
	float	a2;
	float	a;
	float	b;
	float	c;
	float	delta;
	float	t1;
	float	t2;

	x = ft_vectorsub(ray->source, ellipsoid->pos);
	r = ellipsoid->radius1 + ellipsoid->radius2;
	// printf("k = %.2f\n", ellipsoid->distance);
	a1 = 2.0 * ellipsoid->distance * ft_dotproduct(ray->direction, ellipsoid->axis);
	a2 = (r * r) + (2.0 * ellipsoid->distance * ft_dotproduct(x, ellipsoid->axis)) - ellipsoid->distance;
	a = (4.0 * r * r  * ft_dotproduct(ray->direction, ray->direction)) - a1 * a1;
	b = 2.0 * (4.0 * r * r * ft_dotproduct(ray->direction, x)) - (a1 * a2);
	c = (4.0 * r * r * ft_dotproduct(x, x)) - a2 * a2;
	delta = (b * b) - (4.0 * a * c);
	if (delta < 0)
		return (0);
	delta = sqrtf(delta);
	t1 = (-b + delta) / (2 * a);
	t2 = (-b - delta) / (2 * a);
	return (ft_min_ray(t1, t2, tmin));
}

int		paraboloid_intersect(t_object *paraboloid, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	a;
	float	b;
	float	c;
	float	delta;
	float	t;
	float	t1;
	float	t2;

	t = INFINITY;
	x = ft_vectorsub(ray->source, paraboloid->pos);
	// printf("%.2f\n", paraboloid->distance);
	// exit (0);
	a = ft_dotproduct(ray->direction, ray->direction) - powf(ft_dotproduct(ray->direction, paraboloid->axis), 2.0);
	b = 2.0 * (ft_dotproduct(ray->direction, x) - ft_dotproduct(ray->direction, ft_vectormulti(paraboloid->axis, (ft_dotproduct(x, paraboloid->axis) + (2.0 * paraboloid->distance)))));
	c = ft_dotproduct(x, x) - ft_dotproduct(x, ft_vectormulti(paraboloid->axis, (ft_dotproduct(x, paraboloid->axis) + (4.0 * paraboloid->distance))));
	delta = (b * b) - (4.0 * a * c);
	if (delta < 0)
		return (0);
	delta = sqrtf(delta);
	t1 = (-b + delta) / (2 * a);
	t2 = (-b - delta) / (2 * a);
	return (ft_min_ray(t1, t2, tmin));
}

int		triangle_intersect(t_object *triangle, t_ray *ray, float *tmin)
{
	t_vec	x;
	t_vec	ca;
	t_vec	ba;
	t_vec	bc;
	t_vec	ab;
	t_vec	normal;
	t_vec	q;
	t_vec	qa;
	t_vec	qb;
	t_vec	qc;
	float	distance;
	float	dist2plane;
	float	a;
	float	b;
	float	t;
	float	t1;
	float	t2;

	ca = ft_vectorsub(triangle->c, triangle->a);
	ba = ft_vectorsub(triangle->b, triangle->a);
	normal = ft_normalize(ft_crossproduct(ca, ba));
	distance = ft_dotproduct(normal, triangle->a);
	// a = ft_dotproduct(ray->direction, normal);
	// b = ft_dotproduct(normal, ft_vectoradd(ray->source, ft_negative(ft_vectormulti(normal, distance))));
	x = ft_vectorsub(ray->source, triangle->a);
	a = -1.0 * ft_dotproduct(x, normal);
	b = ft_dotproduct(ray->direction, normal);
	dist2plane = a / b;
	q = ft_vectoradd(ray->source, ft_vectormulti(ray->direction, dist2plane));
	ca = ft_vectorsub(triangle->c, triangle->a);
	qa = ft_vectorsub(q, triangle->a);
	t = ft_dotproduct(ft_crossproduct(ca, qa), normal);
	bc = ft_vectorsub(triangle->b, triangle->c);
	qc = ft_vectorsub(q, triangle->c);
	t1 = ft_dotproduct(ft_crossproduct(bc, qc), normal);
	ab = ft_vectorsub(triangle->a, triangle->b);
	qb = ft_vectorsub(q, triangle->b);
	t2 = ft_dotproduct(ft_crossproduct(ab, qb), normal);
	if (t >= 0.0 && t1 >= 0.0 && t2 >= 0.0)
	{
		*tmin = dist2plane;
		return (dist2plane);
	}
	return (0);
}

// int		ltd_pln_intersect(t_object *ltd_pln, t_ray *ray, float *tmin)
// {
// 	t_object	tr1;
// 	t_object	tr2;
// 	t_vec		a;
// 	t_vec		b;
// 	t_vec		c;
// 	int			ret1;
// 	int			ret2;
// 	// t_vec		d;

// 	// d = ft_vector(ltd_pln->corner2.x, ltd_pln->corner1.y, ltd_pln->corner1.z);
// 	a = ft_vector(ltd_pln->corner2.x, ltd_pln->corner2.y, ltd_pln->corner1.z);
// 	b = ft_vector(ltd_pln->corner1.x, ltd_pln->corner2.y, ltd_pln->corner1.z);
// 	c = ft_vector(ltd_pln->corner1.x, ltd_pln->corner2.y, ltd_pln->corner2.z);
// 	tr1.a = a;
// 	tr1.b = b;
// 	tr1.c = c;
// 	tr1.color = ltd_pln->color;
// 	tr1.rot = ltd_pln->rot;
// 	tr1.trans = ltd_pln->trans;
// 	tr2.a = c;
// 	tr2.b = ltd_pln->corner2;
// 	tr2.c = a;
// 	tr2.color = ltd_pln->color;
// 	tr2.rot = ltd_pln->rot;
// 	tr2.trans = ltd_pln->trans;
// 	ret1 = triangle_intersect(&tr1, ray, tmin);
// 	ret2 = triangle_intersect(&tr2, ray, tmin);
// 	return (1);
// }

int		box_intersect(t_object *box, t_ray *ray, float *tmin)
{
	t_vec	t[3];
	int		sign[3];
	t_vec	bounds[2];

	bounds[0] = box->bounds[0];
	bounds[1] = box->bounds[1];
	t[2] = (t_vec){1.0 / ray->direction.x, 1.0 / ray->direction.y, 1.0 / ray->direction.z};
	sign[0] = t[2].x < 0.00001;
	sign[1] = t[2].y < 0.00001;
	sign[2] = t[2].z < 0.00001;
	t[0].x = (bounds[sign[0]].x - ray->source.x) * t[2].x;
	t[1].x = (bounds[1 - sign[0]].x - ray->source.x) * t[2].x;
	t[0].y = (bounds[sign[1]].y - ray->source.y) * t[2].y;
	t[1].y = (bounds[1 - sign[1]].y - ray->source.y) * t[2].y;
	if ((t[0].x > t[1].y) || (t[0].y > t[1].x))
		return (0);
	t[0].x = (t[0].y > t[0].x) ? t[0].y : t[0].x;
	t[1].x = (t[1].y < t[1].x) ? t[1].y : t[1].x;
	t[0].z = (bounds[sign[2]].z - ray->source.z) * t[2].z;
	t[1].z = (bounds[1 - sign[2]].z - ray->source.z) * t[2].z;
	if ((t[0].x > t[1].z) || (t[0].z > t[1].x))
		return (0);
	t[0].x = (t[0].z > t[0].x) ? t[0].z : t[0].x;
	t[1].x = (t[1].z < t[1].x) ? t[1].z : t[1].x;
	return (ft_min_ray(t[0].x, t[1].x, tmin));
}

static int	ft_aux_parallelogram(t_aux_parallelogram *a, t_object *para, t_ray *ray)
{
	if (a->a < 1.0e-6)
		return (0);
	a->q = ft_crossproduct(a->t, a->e01);
	a->b = ft_dotproduct(ray->direction, a->q) * a->invdet;
	if (a->b < 1.0e-6)
		return (0);
	if ((a->a + a->b) > 1.0001f)
	{
		a->e23 = ft_vectorsub(para->d, para->c);
		a->e21 = ft_vectorsub(para->b, para->c);
		a->p2 = ft_crossproduct(ray->direction, a->e21);
		a->det2 = ft_dotproduct(a->e23, a->p2);
		if (fabs(a->det2) < 1.0e-6)
			return (0);
		a->invdet2 = 1.0 / a->det2;
		a->t2 = ft_vectorsub(ray->source, para->c);
		a->a2 = ft_dotproduct(a->t2, a->p2) * a->invdet2;
		if (a->a2 < 1.0e-6)
			return (0);
		a->q2 = ft_crossproduct(a->t2, a->e23);
		a->b2 = ft_dotproduct(ray->direction, a->q2) * a->invdet2;
		if (a->b2 < 1.0e-6)
			return (0);
	}
	return (1);
}

int			parallelogram_intersect(t_object *para, t_ray *ray,float *tmin)
{
	t_aux_parallelogram	a;

	a.e01 = ft_vectorsub(para->b, para->a);
	a.e03 = ft_vectorsub(para->d, para->a);
	a.p = ft_crossproduct(ray->direction, a.e03);
	a.det = ft_dotproduct(a.e01, a.p);
	if (fabs(a.det) >= 1.0e-6)
	{
		a.invdet = 1.0 / a.det;
		a.t = ft_vectorsub(ray->source, para->a);
		a.a = ft_dotproduct(a.t, a.p) * a.invdet;
		if (!ft_aux_parallelogram(&a, para, ray))
			return (0);
		a.t1 = ft_dotproduct(a.e03, a.q) * a.invdet;
		if (a.t1 < *tmin && a.t1 > 1.0e-6)
		{
			*tmin = a.t1;
			return (1);
		}
	}
	return (0);
}

void	ft_cylinder_normal(t_hit *hit, t_ray *ray, int ret)
{
	float	m;

	m = .0;
	if (ret < 0)
		hit->n = ft_vectormulti(hit->object->axis, -1.0);
	else if (ret > 1)
		hit->n = hit->object->axis;
	else
	{
	m = ft_dotproduct(ray->direction, hit->object->axis) * hit->t;
	m += ft_dotproduct(ray->source, hit->object->axis);
	hit->n = ft_vectorsub(hit->p, ft_vectormulti(hit->object->axis, m));
	}
}

t_vec	ft_normal_cone(t_object *cone, t_vec p)
{
	p = ft_vectormulti( p, 1.0 / cone->angle);
	return (ft_vector(p.x, 0.001, p.z));
}

t_vec	ft_box_normal(t_object *box, t_vec h)
{
	t_vec	c;
	t_vec	d;
	t_vec	p;
	float	bias;

	bias = 1.00001;
	c = ft_vectormulti(ft_vectoradd(box->bounds[0], box->bounds[1]), 0.5);
	d = ft_vectormulti(ft_vectorsub(box->bounds[0], box->bounds[1]), 0.5);
	d.x = fabs(d.x) * bias;
	d.y = fabs(d.y) * bias;
	d.z = fabs(d.z) * bias;
	p = ft_vectorsub(h, c);
	return (ft_normalize(ft_vector(p.x / d.x, p.y / d.y, p.z / d.z)));
}

t_vec	ft_parallelogram_normal(t_object *para)
{
	t_vec	n;

	n = ft_crossproduct(ft_vectorsub(para->c, para->a), ft_vectorsub(para->d, para->a));
	return (ft_normalize(n));	
}

void	ft_compute_normals(t_hit *hit, t_ray *ray)
{
	t_vec	x;
	t_vec	cmid;
	t_vec	mr;
	t_vec	ba;
	t_vec	ca;
	// t_vec	ab;
	// t_vec	cb;
	float	m;
	// float	maxm;
	// float	k;
	int		ret;
	
	ret = 0;
	if (hit->object->type == SPHERE)
		hit->n = ft_normalize(ft_vectorsub(hit->p, hit->object->pos));
	else if (hit->object->type == PLANE)
		hit->n = hit->object->normal;
	else if (hit->object->type == CYLINDER)
	{
		if (hit->object->height > 0.0)
		{
			if (!(ret = cylinder_intersect(hit->object, ray, &hit->t)))
				return ;
			ft_cylinder_normal(hit, ray, ret);
		}
		else
		{
			x = ft_vectorsub(ray->source, hit->object->pos);
			m = ft_dotproduct(ray->direction, ft_vectormulti(hit->object->axis, hit->t)) + ft_dotproduct(x, hit->object->axis);
			hit->n = ft_normalize(ft_vectorsub(ft_vectorsub(hit->p, hit->object->pos), ft_vectormulti(hit->object->axis, m)));
		}
	}
	else if (hit->object->type == CONE)
	{
		hit->n = ft_normal_cone(hit->object, hit->p);
		// x = ft_vectorsub(ray->source, hit->object->pos); 
		// m = ft_dotproduct(ray->direction, ft_vectormulti(hit->object->axis, hit->t)) + ft_dotproduct(x, hit->object->axis);
		// k = tanf(deg_to_rad(hit->object->angle) / 2.0);
		// hit->n = ft_normalize(ft_vectorsub(ft_vectorsub(hit->p, hit->object->pos), ft_vectormulti(hit->object->axis, ((1.0 + (k * k)) * m))));
	}
	else if (hit->object->type == ELLIPSOID)
	{
		cmid = ft_vectoradd(hit->object->pos, ft_vectormulti(hit->object->axis, (hit->object->distance / 2.0)));
		mr = ft_vectorsub(hit->p, cmid);
		hit->n = ft_normalize(ft_vectorsub(mr, ft_vectormulti(hit->object->axis, ((1.0f - powf(hit->object->radius2, 2.0) / powf(hit->object->radius1, 2.0)) * ft_dotproduct(mr, hit->object->axis)))));
	}
	else if (hit->object->type == PARABOLOID)
	{
		m = ft_dotproduct(ft_vectorsub(hit->p, hit->object->pos), hit->object->axis);
		hit->n = ft_normalize(ft_vectorsub(ft_vectorsub(hit->p, hit->object->pos), ft_vectormulti(hit->object->axis, (m + hit->object->distance))));
	}
	else if (hit->object->type == TRIANGLE)
	{
		ca = ft_vectorsub(hit->object->c, hit->object->a);
		ba = ft_vectorsub(hit->object->b, hit->object->a);
		hit->n = ft_normalize(ft_crossproduct(ba, ca));
	}
	else if (hit->object->type == BOX)
	{
		hit->n = ft_box_normal(hit->object, hit->p);
	}
	// else if (hit->object->type == LTD_PLANE)
	// {
	// 	hit->n = ltd_pln_normal();
	// }
}

t_vec	ft_light_computing(t_light *light, t_vec light_dir, t_hit *hit, t_ray *ray)
{
	float 		lambert;
	float		ambient_strenght;
	t_vec		ambient;
	t_vec		color;
	float		reflect;
	float		phong_term;
	t_vec		phong_dir;
	// t_vec		blinn_dir;
	// float		temp;
	// float		blinn_term;

	// this is the ambient calculation;

	color = (t_vec){0.0f, 0.0f, 0.0f};
	ambient_strenght = 0.05f;
	ambient = ft_vectormulti(light->color, ambient_strenght);
	color.x = ambient.x * hit->object->color.x;
	color.y = ambient.y * hit->object->color.y;
	color.z = ambient.z * hit->object->color.z;

	// the following is the diffuse calculation;

	//color = (t_vec){0.0f, 0.0f, 0.0f};
	lambert = fmax(0.0f, ft_dotproduct(light_dir, hit->n));
	color = ft_vectoradd(color, ft_vectormulti(hit->object->color, lambert));

	// and this is the calculation of the shininess of the object (specular) using Bui Tuong Phong shading methode;

	reflect = 2.0f * (ft_dotproduct(light_dir, hit->n));
	phong_dir = ft_vectorsub(light_dir, ft_vectormulti(hit->n, reflect));
	phong_term = fmax(ft_dotproduct(phong_dir, ray->direction), 0.0f);
	phong_term = 1.0f * powf(phong_term, 90.0f) * 1.0f;
	color = ft_vectoradd(color, ft_vectormulti(light->color, phong_term));
	
	// and this is the calculation of the shininess of the object (specular) using Jim Blinn shading way;

	// blinn_dir = ft_vectorsub(light_dir, ray->direction);
	// temp = sqrtf(ft_dotproduct(blinn_dir, blinn_dir));
	// if (temp != 0.0f)
	// {
	// 	blinn_dir = ft_vectormulti(blinn_dir, (1.0f / temp));
	// 	blinn_term = fmax(ft_dotproduct(blinn_dir, hit->n), 0.0f);
	// 	blinn_term = 1.0f * powf(blinn_term, 90.0f) * 1.0f;
	// 	color = ft_vectoradd(color, ft_vectormulti(hit->object->color, blinn_term));
	// }
	//printf("Intensity : %.10f\n", light->intensity);
	color.x = color.x * light->color.x * light->intensity;
	color.y = color.y * light->color.y * light->intensity;
	color.z = color.z * light->color.z * light->intensity;
	return (color);
}

int 	ft_shade_object(t_hit *hit, t_light *lights, t_object *lst, t_ray *ray)
{
	t_light 	*light;
	t_vec		color;
	t_vec		light_dir;
	t_vec       dist;
	t_ray		shadow_ray;

	float t;

	color = (t_vec){0.0f, 0.0f, 0.0f};
	shadow_ray.source = hit->p;//ft_vectoradd(hit->p, (t_vec){0.5 * hit->n.x, 0.5 * hit->n.y, 0.5 * hit->n.z});
	light = lights;
	while (light)
	{
		dist = ft_vectorsub(light->pos, hit->p);
		light_dir = ft_normalize(dist);
		// if (ft_dotproduct(hit->n, dist)<=0.0f){
		// 	ft_putendl("Yeah");
		// 	light = light->next;
		// 	continue;
		// 	break;
		// }
		shadow_ray.direction = light_dir;//(t_vec){-light_dir.x, -light_dir.y, -light_dir.z};
		t = ft_magnitude(ft_vectorsub(light->pos, hit->p));
		// if (t <= 0.0f){
		// 	ft_putendl("Yeoh");
		// 	light = light->next;
		// 	continue;
		// 	break;
		// }
		if (!shadow_cast(lst, &shadow_ray, &t))
			color = ft_vectoradd(color, ft_light_computing(light, light_dir, hit, ray));
		light = light->next;
	}
	return (rgb_to_int(clamp_vect(color)));
}

t_vec	ft_rotate_object(t_vec to_rot, t_vec rot, int invert)
{
	if (invert)
	{
		to_rot = z_rotation(to_rot, -rot.z);
		to_rot = y_rotation(to_rot, -rot.y);
		// result = z_rotation(to_rot, -rot.z);
		// result = y_rotation(to_rot, -rot.y);
		// result = x_rotation(to_rot, -rot.x);
		to_rot = x_rotation(to_rot, -rot.x);
	}
	else
	{
		to_rot = x_rotation(to_rot, rot.x);
		to_rot = y_rotation(to_rot, rot.y);
		// result = x_rotation(to_rot, rot.x);
		// result = y_rotation(to_rot, rot.y);
		// result = z_rotation(to_rot, rot.z);
		to_rot = z_rotation(to_rot, rot.z);
	}
	return (to_rot);
}

t_vec	ft_translate_object(t_vec to_trans, t_vec trans, int invert)
{
	t_vec	result;

	if (invert)
	{
		result.x = to_trans.x - trans.x;
		result.y = to_trans.y - trans.y;
		result.z = to_trans.z - trans.z;
	}
	else
	{
		result.x = to_trans.x + trans.x;
		result.y = to_trans.y + trans.y;
		result.z = to_trans.z + trans.z;
	}
	return (result);
}

t_ray	ft_transform_ray(t_object *obj, t_ray *raw, int invert)
{
	t_ray	result;

	result = *raw;
	result.source = ft_rotate_object(raw->source, obj->rot, invert);
	result.source = ft_translate_object(result.source, obj->trans, invert);
	result.direction = ft_rotate_object(raw->direction, obj->rot, invert);
	return (result);
}

int		raycast(t_object *lst, t_ray *raw, t_hit *hit)
{
	t_object	*p;
	float		t;
	t_ray		ray;
	t_ray		save;

	t = INFINITY;
	hit->object = NULL;
	hit->t = INFINITY;
	p = lst;
	while (p)
	{
		ray = ft_transform_ray(p, raw, 1);
		//ray = *raw;
		if (p->type == SPHERE)
		{
			if (sphere_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		else if (p->type == CYLINDER)
		{
			if (cylinder_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		else if (p->type == PLANE)
		{
			if (plane_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		else if (p->type == CONE)
		{
			if (cone_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		else if (p->type == ELLIPSOID)
		{
			if (ellipsoid_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		else if (p->type == PARABOLOID)
		{
			if (paraboloid_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		else if (p->type == TRIANGLE)
		{
			if (triangle_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		else if (p->type == BOX)
		{
			if (box_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		// else if (p->type == LTD_PLANE)
		// {
		// 	if (ltd_pln_intersect(p, &ray, &t))
		// 		if (hit->t > t)
		// 		{
		// 			hit->t = t;
		// 			hit->object = p;
		// 			save = ray;
		// 		}
		// }
		else if (p->type == PARALLELOGRAM)
		{
			if (parallelogram_intersect(p, &ray, &t))
				if (hit->t > t)
				{
					hit->t = t;
					hit->object = p;
					save = ray;
				}
		}
		p = p->next;
	}
	if (hit->object == NULL)
		return (0);
	hit->p = ft_vectoradd(save.source, ft_vectormulti(save.direction, hit->t));
	ft_compute_normals(hit, &save);
	hit->p = ft_translate_object(hit->p, hit->object->trans, 0);
	hit->p = ft_rotate_object(hit->p, hit->object->rot, 0);
	hit->n = ft_rotate_object(hit->n, hit->object->rot, 0);
	hit->n = ft_normalize(hit->n);
	return (1);
}

int 	shadow_cast(t_object *lst, t_ray *ray, float *tmin)
{
	t_object	*p;
	float		t;
	t_ray		ra;

	t = INFINITY;
	p = lst;
	while (p)
	{
		ra = ft_transform_ray(p, ray, 1);
		//ra = *ray;
		if (p->type == SPHERE)
		{
			if (sphere_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		else if (p->type == PLANE)
		{
			if (plane_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		else if (p->type == CYLINDER)
		{
			if (cylinder_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		else if (p->type == CONE)
		{
			if (cone_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		else if (p->type == ELLIPSOID)
		{
			if (ellipsoid_intersect(p, &ra, &t))
			{
				if (t < *tmin)
					return (1);
			}		
		}
		else if (p->type == PARABOLOID)
		{
			if (paraboloid_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		else if (p->type == TRIANGLE)
		{
			if (triangle_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		else if (p->type == BOX)
		{
			if (box_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		// else if (p->type == LTD_PLANE)
		// {
		// 	if (ltd_pln_intersect(p, &ra, &t))
		// 		if (t < *tmin)
		// 			return (1);
		// }
		else if (p->type == PARALLELOGRAM)
		{
			if (parallelogram_intersect(p, &ra, &t))
				if (t < *tmin)
					return (1);
		}
		p = p->next;
	}
	return 0;
}

// void  ft_print_vect(t_vec v){
// 	printf("x : %.2f y : %.2f z : %.2f\n", v.x, v.y, v.z);
// }


void update(t_mx *mx)
{
	t_ray 	ray;
	t_hit	hit;
	int		y;
	int		x;
 
    x = -1;
	while (++x < WIN_W)
	{
		y = -1;
		while (++y < WIN_H)
		{

			ray = camera_ray(&mx->cam, x, y);
			hit.t = INFINITY;
			if (raycast(mx->objects, &ray, &hit))
			{
				mx->rt[x + (WIN_H - 1 - y) * WIN_W] = ft_shade_object(&hit, mx->lights, mx->objects, &ray);
				// mx->rt[x + (WIN_H - 1 - y) * WIN_W] = 0xff0000;
			}
		}
	}
}

int     main(int ac, char **av)
{
    t_mx    v;

    if (ac == 2)
    {
        if (av[1])
        {
			ft_bzero(&v, sizeof(t_mx));
			//ft_memset(&v, 0, sizeof(t_mx));
            if (!ft_open(av[1], &v)){
                ft_putstr("error! please try a valid configuration file.\n");
				exit(0);
            }
            else
				run(&v);
        }
    }
    else
        ft_usage();
    return (0);
}