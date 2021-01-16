/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   objects.c                                          :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/10/14 18:10:13 by ahkhilad          #+#    #+#             */
/*   Updated: 2021/01/16 19:05:00 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

t_object        *ft_object_new(t_object o)
{
    t_object    *new;

    if ((new = ft_memalloc(sizeof(t_object))))
    {
        new->type = o.type;
        new->pos = o.pos;
        new->a = o.a;
        new->b = o.b;
        new->c = o.c;
        new->d = o.d;
        // new->corner1 = o.corner1;
        // new->corner2 = o.corner2;
        new->bounds[0] = o.bounds[0];
        new->bounds[1] = o.bounds[1];
        new->height = o.height;
        new->trans = o.trans;
        new->rot = o.rot;
        new->radius = o.radius;
        new->radius1 = o.radius1;
        new->radius2 = o.radius2;
        new->distance = o.distance;
        new->angle = o.angle;
        new->axis = o.axis;
        new->normal = o.normal;
        new->color = o.color;
        new->next = NULL;
    }
    return (new);
}

void        ft_object_push(t_object **lst, t_object *new)
{
    t_object    *p;

    if (lst == NULL || new == NULL)
    return;
    else if (*lst == NULL)
    {
        *lst = new;
        return ;
    }
    p = *lst;
    while (p->next != NULL)
        p = p->next;
    p->next = new;
}

void    ft_object_clear(t_object **lst)
{
    t_object    *p;

    if (lst == NULL)
        return ;
    while (*lst)
    {
        p = *lst;
        if (p->next)
            *lst = p->next;
        free(p);
    }
}