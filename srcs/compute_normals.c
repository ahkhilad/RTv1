/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   compute_normals.c                                  :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/12/27 18:48:06 by babdelka          #+#    #+#             */
/*   Updated: 2021/01/23 23:41:27 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

void	ft_cylinder_normal(t_hit *hit, t_ray *ray, int ret)
{
	float	m;

	m = 0.0;
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

t_vec		ft_torus_normal(t_object *torus, t_vec h)
{
	t_vec	p;
	float	x;
	float	y;
	float	z;

	x = h.x;
	y = h.y;
	z = h.z;
	p.x = 4.0 * x * (x * x + y * y + z * z - torus->radius2 * torus->radius2
			- torus->radius1 * torus->radius1);
	p.y = 4.0 * y * (x * x + y * y + z * z - torus->radius2 * torus->radius2
			- torus->radius1 * torus->radius1 + 2.0 * torus->radius1 * torus->radius1);
	p.z = 4.0 * z * (x * x + y * y + z * z - torus->radius2 * torus->radius2
			- torus->radius1 * torus->radius1);
	// return (ft_normalize(p));
	return (p);
}

void		ft_compute_normals(t_hit *hit, t_ray *ray)
{
	t_vec	x;
	float	m;
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
		/*x = ft_vectorsub(ray->source, hit->object->pos);
		 ft_dotproduct(ray->direction, ft_vectormulti(hit->object->axis,\
		->t)) + ft_dotproduct(x, hit->object->axis);
		 tanf(deg_to_rad(hit->object->angle) / 2.0);
		->n = ft_normalize(ft_vectorsub(ft_vectorsub(hit->p, hit->object->\
		), ft_vectormulti(hit->object->axis, ((1.0 + (k * k)) * m))));*/
	}
	else if (hit->object->type == BOX)
	{
		hit->n = ft_box_normal(hit->object, hit->p);
	}
	else if (hit->object->type == PARALLELOGRAM)
	{
		hit->n = ft_parallelogram_normal(hit->object);
	}
	else if (hit->object->type == TORUS)
	{
		hit->n = ft_torus_normal(hit->object, hit->p);
	}
}
