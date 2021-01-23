/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   tokenfill.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/12/26 18:54:28 by babdelka          #+#    #+#             */
/*   Updated: 2021/01/23 23:21:54 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

void	tokenfill_sphere(char **token, t_object *object)
{
	if (token[1])
		object->pos = string_to_vect(token[1]);
	if (token[2])
		object->radius = ft_atof(token[2]);
	if (token[3])
		object->rot = string_to_vect(token[3]);
	if (token[4])
		object->trans = string_to_vect(token[4]);
	if (token[5])
		object->color = vect_from_hexa(ft_special_atoi_base(token[5]));
}

void	tokenfill_plane(char **token, t_object *object)
{
	if (token[1])
		object->pos = string_to_vect(token[1]);
	if (token[2])
		object->normal = string_to_vect(token[2]);
	if (token[3])
		object->rot = string_to_vect(token[3]);
	if (token[4])
		object->trans = string_to_vect(token[4]);
	if (token[5])
		object->color = vect_from_hexa(ft_special_atoi_base(token[5]));
}

void	tokenfill_cone(char **token, t_object *object)
{
	if (token[1])
		object->pos = string_to_vect(token[1]);
	if (token[2])
		object->angle = ft_atof(token[2]);
	if (token[3])
		object->height = ft_atof(token[3]);
	if (token[4])
		object->axis = string_to_vect(token[4]);
	if (token[5])
		object->rot = string_to_vect(token[5]);
	if (token[6])
		object->trans = string_to_vect(token[6]);
	if (token[7])
		object->color = vect_from_hexa(ft_special_atoi_base(token[7]));
}

void	tokenfill_cylinder(char **token, t_object *object)
{
	if (token[1])
		object->pos = string_to_vect(token[1]);
	if (token[2])
		object->radius = ft_atof(token[2]);
	if (token[3])
		object->height = ft_atof(token[3]);
	if (token[4])
		object->axis = string_to_vect(token[4]);
	if (token[5])
		object->rot = string_to_vect(token[5]);
	if (token[6])
		object->trans = string_to_vect(token[6]);
	if (token[7])
		object->color = vect_from_hexa(ft_special_atoi_base(token[7]));
}

void	tokenfill_box(char **token, t_object *object)
{
	if (token[1])
        object->bounds[0] = string_to_vect(token[1]);
    if (token[2])
        object->bounds[1] = string_to_vect(token[2]);
    if (token[3])
        object->rot = string_to_vect(token[3]);
    if (token[4])
        object->trans = string_to_vect(token[4]);
    if (token[5])
        object->color = vect_from_hexa(ft_special_atoi_base(token[5]));
}

void	tokenfill_parallelogram(char **token, t_object *object)
{
    if (token[1])
        object->a = string_to_vect(token[1]);
    if (token[2])
        object->b = string_to_vect(token[2]);
    if (token[3])
        object->c = string_to_vect(token[3]);
    if (token[4])
        object->d = string_to_vect(token[4]);
    if (token[5])
        object->rot = string_to_vect(token[5]);
    if (token[6])
        object->trans = string_to_vect(token[6]);
    if (token[7])
        object->color = vect_from_hexa(ft_special_atoi_base(token[7]));
}

void	tokenfill_torus(char **token, t_object *object)
{
	if (token[1])
		object->pos = string_to_vect(token[1]);
	if (token[2])
		object->radius1 = ft_atof(token[2]);
	if (token[3])
		object->radius2 = ft_atof(token[3]);
	if (token[4])
		object->rot = string_to_vect(token[4]);
	if (token[5])
		object->trans = string_to_vect(token[5]);
	if (token[6])
		object->color = vect_from_hexa(ft_special_atoi_base(token[6]));
}
