/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   intersection.c                                     :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/12/26 14:28:35 by babdelka          #+#    #+#             */
/*   Updated: 2021/01/23 23:31:11 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

static int	ft_min_ray(float t1, float t2, float *t)
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

int			sphere_intersect(t_object *sphere, t_ray *ray, float *tmin)
{
	t_vec x;

	x = ft_vectorsub(ray->source, sphere->pos);
	sphere->sp.a = ft_dotproduct(ray->direction, ray->direction);
	sphere->sp.b = 2.0 * ft_dotproduct(ray->direction, x);
	sphere->sp.c = ft_dotproduct(x, x) - sphere->radius * sphere->radius;
	sphere->sp.delta = (sphere->sp.b * sphere->sp.b) -\
	(4.0 * sphere->sp.a * sphere->sp.c);
	if (sphere->sp.delta < 0)
		return (0);
	sphere->sp.delta = sqrtf(sphere->sp.delta);
	sphere->sp.t1 = (-sphere->sp.b + sphere->sp.delta) / (2 * sphere->sp.a);
	sphere->sp.t2 = (-sphere->sp.b - sphere->sp.delta) / (2 * sphere->sp.a);
	return (ft_min_ray(sphere->sp.t1, sphere->sp.t2, tmin));
}

int			plane_intersect(t_object *plane, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	a;
	float	b;
	float	t;

	x = ft_vectorsub(ray->source, plane->pos);
	a = -1.0 * ft_dotproduct(x, plane->normal);
	b = ft_dotproduct(ray->direction, plane->normal);
	if (fabs(b) <= 1e-6)
		return (0);
	t = a / b;
	if (t > 1e-2 && t < *tmin)
	{
		*tmin = t;
		return (1);
	}
	return (0);
}

static int	ft_cylinder_cap(t_ray *ray, float *t, t_vec pos, t_vec n)
{
	float	ddn;
	float	t1;
	t_vec	dist;

	ddn = ft_dotproduct(ray->direction, n);
	if (ddn <= 1.0e-6 && ddn >= -1.0e-6)
		return (0);
	dist = ft_vectorsub(pos, ray->source);
	// dist = ft_vectorsub(ray->source, pos);
	t1 = (ft_dotproduct(dist, n)) / ddn;
	//if (t1 < *t && t1 > 1e-2)
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
	if (ft_cylinder_cap(ray, t, ft_vectormulti(cylinder->axis, -cylinder->height), ft_vectormulti(cylinder->axis, -1.0)))
		return (-1);
	return (0);
}

int			cylinder_intersect(t_object *cylinder, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	m0;
	float	m1;

	x = ft_vectorsub(ray->source, cylinder->pos);
	cylinder->cy.a = ft_dotproduct(ray->direction, ray->direction) -\
	pow(ft_dotproduct(ray->direction, cylinder->axis), 2.0);
	cylinder->cy.b = 2.0 * (ft_dotproduct(ray->direction, x) -\
	(ft_dotproduct(ray->direction, cylinder->axis) *\
	ft_dotproduct(x, cylinder->axis)));
	cylinder->cy.c = ft_dotproduct(x, x) -\
	pow(ft_dotproduct(x, cylinder->axis), 2.0) -\
	(cylinder->radius * cylinder->radius);
	cylinder->cy.delta = (cylinder->cy.b * cylinder->cy.b) -\
	(4.0 * cylinder->cy.a * cylinder->cy.c);
	if (cylinder->cy.delta < 0)
		return (0);
	cylinder->cy.delta = sqrt(cylinder->cy.delta);
	cylinder->cy.t1 = (-cylinder->cy.b + cylinder->cy.delta) /\
	(2 * cylinder->cy.a);
	cylinder->cy.t2 = (-cylinder->cy.b - cylinder->cy.delta) /\
	(2 * cylinder->cy.a);
	if (cylinder->height <= 0)
		return (ft_min_ray(cylinder->cy.t1, cylinder->cy.t2, tmin));
	if (cylinder->cy.t1 > cylinder->cy.t2)
	{
		m0 = cylinder->cy.t1;
		cylinder->cy.t1 = cylinder->cy.t2;
		cylinder->cy.t2 = m0;
	}
	m0 = ft_dotproduct(ray->direction, cylinder->axis) * cylinder->cy.t1;
	m0 += ft_dotproduct(x, cylinder->axis);
	m1 = ft_dotproduct(ray->direction, cylinder->axis) * cylinder->cy.t2;
	m1 += ft_dotproduct(x, cylinder->axis);
    //printf("m0 : %.2f\t  m1 : %.2f\n", m0, m1);
	if (m0 < -cylinder->height)
		return (ft_cylinder_cap1(cylinder, ray, tmin, m1));
	else if (m0 >= -cylinder->height && m0 <= cylinder->height)
	{
		if (cylinder->cy.t1 < 1e-6 || cylinder->cy.t1 > *tmin)
			return (0);
		*tmin = cylinder->cy.t1;
		return (1);
	}
	else if (m0 > cylinder->height)
	{
		if (m1 > cylinder->height)
			return (0);
		if (ft_cylinder_cap(ray, tmin, ft_vectormulti(cylinder->axis, cylinder->height), cylinder->axis))
			return (2);
	}
	return (0);
}

float			ft_get_cone_limits(float t, t_ray *ray, t_vec va, t_vec v)
{
	t_vec	q;
	float	m1;
	float	m2;

	if (t < 0)
		return (INFINITY);
	q = ft_vectoradd(ray->source, ft_vectormulti(ray->direction, t));
	m1 = ft_dotproduct(va, q);
	m2 = ft_dotproduct(va, ft_vectorsub(q, v));
	if (m1 < 0.0 && m2 > 0.0)
		return (t);
	return (INFINITY);
}

int			cone_intersect(t_object *cone, t_ray *ray, float *tmin)
{
	t_vec	x;
	t_vec	va;
	t_vec	a;
	t_vec	b;
	float	m;
	float	anglesin;
	float	anglecos;
	// float	s1;
	// float	s2;

	/* x = ft_vectorsub(ray->source, cone->pos);
	va = ft_normalize(ft_vectorsub(cone->pos, cone->axis));
	cone->cn.k = tanf(deg_to_rad(cone->angle) / 2.0);
	cone->cn.a = ft_dotproduct(ray->direction, ray->direction) -\
	(1.0 + (cone->cn.k * cone->cn.k)) *\
	pow(ft_dotproduct(ray->direction, cone->axis), 2.0);
	cone->cn.b = 2.0 * (ft_dotproduct(ray->direction, x) -\
	((1.0 + (cone->cn.k * cone->cn.k)) *\
	ft_dotproduct(ray->direction, cone->axis) * ft_dotproduct(x, cone->axis)));
	cone->cn.c = ft_dotproduct(x, x) -\
	(1.0 + (cone->cn.k * cone->cn.k)) * pow(ft_dotproduct(x, cone->axis), 2.0);*/
	// printf("%.3f\n", cone->height);
	// printf("%.3f\n", cone->angle);
	// exit(0);
	anglecos = pow(cos(cone->angle), 2.0);
	anglesin = pow(sin(cone->angle), 2.0);
	va = ft_normalize(ft_vectorsub(cone->pos, cone->axis));
	x = ft_vectorsub(ray->source, cone->pos);
	a = ft_vectorsub(ray->direction, ft_vectormulti(va, ft_dotproduct(ray->direction, va)));
	b = ft_vectorsub(x, ft_vectormulti(va, ft_dotproduct(x, va)));
	cone->cn.a = anglecos * ft_dotproduct(a, a) -
		anglesin * ft_dotproduct(ray->direction, va) * ft_dotproduct(ray->direction, va);
	cone->cn.b  = 2.0 * anglecos * ft_dotproduct(a, b) -
		2.0 * anglesin * ft_dotproduct(ray->direction, va) * ft_dotproduct(x, va);
	cone->cn.c = anglecos * ft_dotproduct(b, b) -
		anglesin * ft_dotproduct(x, va) * ft_dotproduct(x, va);
	cone->cn.delta = cone->cn.b * cone->cn.b - 4.0 * cone->cn.a * cone->cn.c;
	if (cone->cn.delta < 0.00000001)
		return (0);
	cone->cn.delta = sqrt(cone->cn.delta);
	if (cone->height <= 0)
	{
		return (ft_min_ray((-cone->cn.b + cone->cn.delta) / (2.0 * cone->cn.a),
			(-cone->cn.b - cone->cn.delta) / (2.0 * cone->cn.a), tmin));
	}
	cone->cn.t1 = ft_get_cone_limits((-cone->cn.b + cone->cn.delta) / (2.0 * cone->cn.a),
			ray, va, cone->axis);
	cone->cn.t2 = ft_get_cone_limits((-cone->cn.b - cone->cn.delta) / (2.0 * cone->cn.a),
			ray, va, cone->axis);
	if (!ft_min_ray(cone->cn.t1, cone->cn.t2, tmin))
		return (0);
	m = ft_dotproduct(ray->direction, cone->axis) * (*tmin);
	m += ft_dotproduct(ray->source, cone->axis);
	if (m > cone->height / 2.0 - 5.0e-1)
		return (2);
	return (1);
}

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

static int	ft_parallelogram(t_parallelogram *a, t_object *para, t_ray *ray)
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
	t_parallelogram	a;

	a.e01 = ft_vectorsub(para->b, para->a);
	a.e03 = ft_vectorsub(para->d, para->a);
	a.p = ft_crossproduct(ray->direction, a.e03);
	a.det = ft_dotproduct(a.e01, a.p);
	if (fabs(a.det) >= 1.0e-6)
	{
		a.invdet = 1.0 / a.det;
		a.t = ft_vectorsub(ray->source, para->a);
		a.a = ft_dotproduct(a.t, a.p) * a.invdet;
		if (!ft_parallelogram(&a, para, ray))
			return (0);
		a.t1 = ft_dotproduct(a.e03, a.q) * a.invdet;
		if (a.t1 < *tmin && a.t1 > 1.0e-2)
		// if (a.t1 >= 0.0)
		{
			*tmin = a.t1;
			return (1);
		}
	}
	return (0);
}

/////////////////////////////////////////////////////////////// 					TORUS						///////////////////////////////////////////////////////////////

static int	ft_min_ray_torus(double s[4], float *t, int n)
{
	double	min;
	int		i;

	i = -1;
	while (++i < n)
	{
		if (s[i] > 0.0001)
		{
			min = s[i];
			break ;
		}
	}
	if (i == n)
		return (0);
	while (++i < n)
	{
		if (s[i] > 0.0001 && s[i] < min)
			min = s[i];
	}
	if (min < *t && min > 0.0001)
	{
		*t = min;
		return (n);
	}
	return (0);
}

static int	ft_case_one_cubic(double s[2], double *q)
{
	int		num;
	double	u;

	if (*q > -1e-9 && *q < 1e-9)
	{
		s[0] = 0;
		num = 1;
	}
	else
	{
		u = ((-(*q)) > 0.0 ? (pow((double)(-(*q)), 1.0 / 3.0)) : ((-(*q)) < 0.0 ? -(pow((double)-(-(*q)), 1.0 / 3.0)) : 0.0));
		s[0] = 2.0 * u;
		s[1] = -u;
		num = 2;
	}
	return (num);
}

static int	ft_case_two_cubic(double s[3], double *q, double *p, double *cb_p)
{
	double	phi;
	double	t;

	phi = 1.0 / 3.0 * acos(-(*q) / sqrt(-(*cb_p)));
	t = 2.0 * sqrt(-(*p));
	s[0] = t * cos(phi);
	s[1] = -t * cos(phi + M_PI / 3.0);
	s[2] = -t * cos(phi - M_PI / 3.0);
	return (3);
}

static int	ft_case_three_cubic(double s[3], double *q, double *d)
{
	double	sqrt_d;
	double	u;
	double	v;

	sqrt_d = sqrt(*d);
	u = ((sqrt_d - (*q)) > 0.0 ? (pow((double)(sqrt_d - (*q)), 1.0 / 3.0)) : ((sqrt_d - (*q)) < 0.0 ? -(pow((double)-(sqrt_d - (*q)), 1.0 / 3.0)) : 0.0));
	v = ((sqrt_d + (*q)) > 0.0 ? (pow((double)(sqrt_d + (*q)), 1.0 / 3.0)) : ((sqrt_d + (*q)) < 0.0 ? -(pow((double)-(sqrt_d + (*q)), 1.0 / 3.0)) : 0.0));
	s[0] = u + v;
	return (1);
}

int			ft_solve_cubic(double w[4], double s[3])
{
	int		i;
	int		num;
	double	sub;
	t_cubic	c;

	c.a = w[2] / w[3];
	c.b = w[1] / w[3];
	c.c = w[0] / w[3];
	c.sq_a = c.a * c.a;
	c.p = 1.0 / 3.0 * (-1.0 / 3.0 * c.sq_a + c.b);
	c.q = 1.0 / 2.0 * (2.0 / 27.0 * c.a * c.sq_a - 1.0 / 3.0 * c.a * c.b + c.c);
	c.cb_p = c.p * c.p * c.p;
	c.d = c.q * c.q + c.cb_p;
	if (c.d > -1e-9 && c.d < 1e-9)
		num = ft_case_one_cubic(s, &c.q);
	else if (c.d < 0.0f)
		num = ft_case_two_cubic(s, &c.q, &c.p, &c.cb_p);
	else
		num = ft_case_three_cubic(s, &c.q, &c.d);
	sub = 1.0 / 3.0 * c.a;
	i = -1;
	while (++i < num)
		s[i] -= sub;
	return (num);
}

int		ft_solve_quadric(double c[3], double s[2])
{
	double	p;
	double	q;
	double	d;

	p = c[1] / (2.0 * c[2]);
	q = c[0] / c[2];
	d = p * p - q;
	if (d > -1e-9 && d < 1e-9)
	{
		s[0] = -p;
		return (1);
	}
	else if (d < 0)
		return (0);
	else
	{
		d = sqrt(d);
		s[0] = d - p;
		s[1] = -d - p;
		return (2);
	}
}

static int	ft_case_one_quartic(double coeffs[4], double s[4], double *q, double *p)
{
	int		num;

	coeffs[0] = *q;
	coeffs[1] = *p;
	coeffs[2] = 0.0;
	coeffs[3] = 1.0;
	num = ft_solve_cubic(coeffs, s);
	s[num++] = 0.0;
	return (num);
}

static int	ft_case_norm_quartic(t_quartic *q)
{
	if (q->u > -1e-9 && q->u < 1e-9)
		q->u = 0;
	else if (q->u > 0.0f)
		q->u = sqrt(q->u);
	else
		return (0);
	if (q->v > -1e-9 && q->v < 1e-9)
		q->v = 0.0;
	else if (q->v > 0.0f)
		q->v = sqrt(q->v);
	else
		return (0);
	return (1);
}

static int	ft_case_two_quartic(double coeffs[4], double s[4], t_quartic *q)
{
	int		num;

	coeffs[0] = 1.0 / 2.0 * q->r * q->p - 1.0 / 8.0 * q->q * q->q;
	coeffs[1] = -q->r;
	coeffs[2] = -1.0 / 2.0 * q->p;
	coeffs[3] = 1.0;
	(void)ft_solve_cubic(coeffs, s);
	q->z = s[0];
	q->u = q->z * q->z - q->r;
	q->v = 2.0 * q->z - q->p;
	if (!ft_case_norm_quartic(q))
		return (0);
	coeffs[0] = q->z - q->u;
	coeffs[1] = q->q < 0.0f ? -q->v : q->v;
	coeffs[2] = 1;
	num = ft_solve_quadric(coeffs, s);
	coeffs[0] = q->z + q->u;
	coeffs[1] = q->q < 0.0f ? q->v : -q->v;
	coeffs[2] = 1;
	num += ft_solve_quadric(coeffs, s + num);
	return (num);
}

int			ft_solve_quartic(double w[5], double s[4])
{
	double		coeffs[4];
	t_quartic	q;
	int			i;
	int			num;

	q.a = w[3] / w[4];
	q.b = w[2] / w[4];
	q.c = w[1] / w[4];
	q.d = w[0] / w[4];
	q.sq_a = q.a * q.a;
	q.p = -3.0 / 8.0 * q.sq_a + q.b;
	q.q = 1.0 / 8.0 * q.sq_a * q.a - 1.0 / 2.0 * q.a * q.b + q.c;
	q.r = -3.0 / 256.0 * q.sq_a * q.sq_a + 1.0 / 16.0 * q.sq_a * q.b -
		1.0 / 4.0 * q.a * q.c + q.d;
	if (q.r > -1e-9 && q.r < 1e-9)
		num = ft_case_one_quartic(coeffs, s, &q.q, &q.p);
	else
		num = ft_case_two_quartic(coeffs, s, &q);
	q.sub = 1.0 / 4.0 * q.a;
	i = -1;
	while (++i < num)
		s[i] -= q.sub;
	return (num);
}

int			torus_intersect(t_object *torus, t_ray *ray, float *tmin)
{
	t_vec	dist;
	int		num;
	double	c[5];
	double	s[4];

	c[4] = pow(ft_dotproduct(ray->direction, ray->direction), 2);
	dist = ray->source;
	c[3] = 4.0 * ft_dotproduct(ray->direction, ray->direction) * ft_dotproduct(dist, ray->direction);
	c[2] = 2.0 * ft_dotproduct(ray->direction, ray->direction) * (ft_dotproduct(dist, dist) -
			(torus->radius1 * torus->radius1 + torus->radius2 * torus->radius2)) + 4.0 *
		pow(ft_dotproduct(dist, ray->direction), 2) +
		4.0 * torus->radius1 * torus->radius1 * ray->direction.y * ray->direction.y;
	c[1] = 4.0 * (ft_dotproduct(dist, dist) - (torus->radius1 * torus->radius1 +
				torus->radius2 * torus->radius2)) * ft_dotproduct(dist, ray->direction) +
		8.0 * (torus->radius1 * torus->radius1) * dist.y * ray->direction.y;
	c[0] = pow(ft_dotproduct(dist, dist) - (torus->radius1 * torus->radius1 +
				torus->radius2 * torus->radius2), 2) - 4.0 * (torus->radius1 * torus->radius1) *
		(torus->radius2 * torus->radius2 - dist.y * dist.y);
	num = ft_solve_quartic(c, s);
	if (num == 0)
		return (0);
	else
		return (ft_min_ray_torus(s, tmin, num));
}
