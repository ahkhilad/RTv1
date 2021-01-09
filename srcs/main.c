/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/01/04 21:42:11 by ahkhilad          #+#    #+#             */
/*   Updated: 2021/01/09 19:20:04 by ahkhilad         ###   ########.fr       */
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

int		cone_intersect(t_object *cone, t_ray *ray, float *tmin)
{
	t_vec	x;
	float	k;
	float	a;
	float	b;
	float	c;
	float	delta;
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
}

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

void	ft_compute_normals(t_hit *hit, t_ray *ray)
{
	t_vec	x;
	t_vec	cmid;
	t_vec	mr;
	t_vec	ba;
	t_vec	ca;
	float	m;
	// float	maxm;
	float	k;
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
		x = ft_vectorsub(ray->source, hit->object->pos); 
		m = ft_dotproduct(ray->direction, ft_vectormulti(hit->object->axis, hit->t)) + ft_dotproduct(x, hit->object->axis);
		k = tanf(deg_to_rad(hit->object->angle) / 2.0);
		hit->n = ft_normalize(ft_vectorsub(ft_vectorsub(hit->p, hit->object->pos), ft_vectormulti(hit->object->axis, ((1.0 + (k * k)) * m))));
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
				mx->rt[x + (WIN_H - 1 - y) * WIN_W] = ft_shade_object(&hit, mx->lights, mx->objects, &ray);// 0xFF00FF;
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