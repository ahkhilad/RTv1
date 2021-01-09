/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   parse_objects.c                                    :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/10/14 18:11:31 by ahkhilad          #+#    #+#             */
/*   Updated: 2021/01/08 16:44:58 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

int    ft_parse_sphere(t_mx *v, char **token)
{
    t_object    object;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 6 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_vectors(token[3])
        && ft_check_vectors(token[4]) && ft_check_positive_floats(token[2])
        && ft_check_hex_color(token[5]))
        {
            ft_bzero(&object, sizeof(t_object));
            object.type = SPHERE;
            if (token[1])
                object.pos = string_to_vect(token[1]);
            if (token[2])
                object.radius = ft_atof(token[2]);
            if (token[3])
                object.rot = string_to_vect(token[3]);
            if (token[4])
                object.trans = string_to_vect(token[4]);
            if (token[5])
                object.color = vect_from_hexa(ft_special_atoi_base(token[5]));
            ft_object_push(&v->objects, ft_object_new(object));
        }
        else
            return (0);   
    }
    else
        return (0);
    return (1);
}

int    ft_parse_plane(t_mx *v, char **token)
{
    t_object    object;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 6 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_directional_vec(token[2])
        && ft_check_vectors(token[3]) && ft_check_vectors(token[4])
        && ft_check_hex_color(token[5]))
        {
            ft_bzero(&object, sizeof(t_object));
            object.type = PLANE;
            if (token[1])
                object.pos = string_to_vect(token[1]);
            if (token[2])
                object.normal = string_to_vect(token[2]);
            if (token[3])
                object.rot = string_to_vect(token[3]);
            if (token[4])
                object.trans = string_to_vect(token[4]);
            if (token[5])
                object.color = vect_from_hexa(ft_special_atoi_base(token[5]));
            ft_object_push(&v->objects, ft_object_new(object));
        }
        else
            return (0);
    }
    else
        return (0);
    return (1);
}

int    ft_parse_cone(t_mx *v, char **token)
{
    t_object    object;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 7 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_directional_vec(token[3])
        && ft_check_vectors(token[4]) && ft_check_vectors(token[5])
        && ft_check_floats(token[2]) && ft_check_hex_color(token[6]))
        {
            ft_bzero(&object, sizeof(t_object));
            object.type = CONE;
            if (token[1])
                object.pos = string_to_vect(token[1]);
            if (token[2])
                object.angle = ft_atof(token[2]);
            if (token[3])
                object.axis = string_to_vect(token[3]);
            if (token[4])
                object.rot = string_to_vect(token[4]);
            if (token[5])
                object.trans = string_to_vect(token[5]);
            if (token[6])
                object.color = vect_from_hexa(ft_special_atoi_base(token[6]));
            ft_object_push(&v->objects, ft_object_new(object));
        }
        else
            return (0);
    }
    else
        return (0);
    return (1);
}

int    ft_parse_cylinder(t_mx *v, char **token)
{
    t_object    object;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 8 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_directional_vec(token[4])
        && ft_check_vectors(token[5]) && ft_check_vectors(token[6])
        && ft_check_positive_floats(token[2]) && ft_check_floats(token[3]) 
        && ft_check_hex_color(token[7]))
        {
            ft_bzero(&object, sizeof(t_object));
            object.type = CYLINDER;
            if (token[1])
                object.pos = string_to_vect(token[1]);
            if (token[2])
                object.radius = ft_atof(token[2]);
            if (token[3])
                object.height = ft_atof(token[3]);
            if (token[4])
                object.axis = string_to_vect(token[4]);
            if (token[5])
                object.rot = string_to_vect(token[5]);
            if (token[6])
                object.trans = string_to_vect(token[6]);
            if (token[7])
                object.color = vect_from_hexa(ft_special_atoi_base(token[7]));
            ft_object_push(&v->objects, ft_object_new(object));
        }
        else
            return (0);
    }
    else
        return (0);
    return (1);
}

int    ft_parse_ellipsoid(t_mx *v, char **token)
{
    t_object    object;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 9 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_positive_floats(token[2])
        && ft_check_positive_floats(token[3]) && ft_check_positive_floats(token[4])
        && ft_check_directional_vec(token[5]) && ft_check_vectors(token[6]) 
        && ft_check_vectors(token[7]) && ft_check_hex_color(token[8]))
        {
            ft_bzero(&object, sizeof(t_object));
            object.type = ELLIPSOID;
            if (token[1])
                object.pos = string_to_vect(token[1]);
            if (token[2])
                object.distance = ft_atof(token[2]);
            if (token[3])
                object.radius1 = ft_atof(token[3]);
            if (token[4])
                object.radius2 = ft_atof(token[4]);
            if (token[5])
                object.axis = string_to_vect(token[5]);
            if (token[6])
                object.rot = string_to_vect(token[6]);
            if (token[7])
                object.trans = string_to_vect(token[7]);
            if (token[8])
                object.color = vect_from_hexa(ft_special_atoi_base(token[8]));
            ft_object_push(&v->objects, ft_object_new(object));
        }
        else
            return (0);
    }
    else
        return (0);
    return (1);
}

int    ft_parse_paraboloid(t_mx *v, char **token)
{
    t_object    object;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 7 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_directional_vec(token[3])
        && ft_check_vectors(token[4]) && ft_check_vectors(token[5])
        && ft_check_floats(token[2]) && ft_check_hex_color(token[6]))
        {
            ft_bzero(&object, sizeof(t_object));
            object.type = PARABOLOID;
            if (token[1])
                object.pos = string_to_vect(token[1]);
            if (token[2])
                object.distance = ft_atof(token[2]);
            if (token[3])
                object.axis = string_to_vect(token[3]);
            if (token[4])
                object.rot = string_to_vect(token[4]);
            if (token[5])
                object.trans = string_to_vect(token[5]);
            if (token[6])
                object.color = vect_from_hexa(ft_special_atoi_base(token[6]));
            ft_object_push(&v->objects, ft_object_new(object));
        }
        else
            return (0);
    }
    else
        return (0);
    return (1);
}

int    ft_parse_triangle(t_mx *v, char **token)
{
    t_object    object;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 7 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_vectors(token[2])
        && ft_check_vectors(token[3]) && ft_check_vectors(token[4])
        && ft_check_vectors(token[5]) && ft_check_hex_color(token[6]))
        {
            ft_bzero(&object, sizeof(t_object));
            object.type = TRIANGLE;
            if (token[1])
                object.a = string_to_vect(token[1]);
            if (token[2])
                object.b = string_to_vect(token[2]);
            if (token[3])
                object.c = string_to_vect(token[3]);
            if (token[4])
                object.rot = string_to_vect(token[4]);
            if (token[5])
                object.trans = string_to_vect(token[5]);
            if (token[6])
                object.color = vect_from_hexa(ft_special_atoi_base(token[6]));
            ft_object_push(&v->objects, ft_object_new(object));
        }
        else
            return (0);
    }
    else
        return (0);
    return (1);
}

int    ft_parse_light(t_mx *v, char **token)
{
    t_light    light;
    int         len;

    len = ft_strsplit_len(token);
    if (len == 4 && token)
    {
        if (ft_check_vectors(token[1]) && ft_check_positive_floats(token[2])
        && ft_check_hex_color(token[3]))
        {
            ft_bzero(&light, sizeof(t_light));
            if (token[1])
                light.pos = string_to_vect(token[1]);
            if (token[2])
            {
                light.intensity = ft_atof(token[2]);
            }
            if (token[3])
                light.color = vect_from_hexa(ft_special_atoi_base(token[3]));
            ft_light_push(&v->lights, ft_light_new(light.pos, light.color, light.intensity));
        }
        else
            return (0);
    }
    else
        return (0);
    return (1);
}