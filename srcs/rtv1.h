/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   rtv1.h                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <marvin@42.fr>                    +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/11/29 00:31:07 by ahkhilad          #+#    #+#             */
/*   Updated: 2019/11/29 00:31:12 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef RTV1_H
# define RTV1_H

# define BUFF_SIZE 100000

# define WIN_W 1280
# define WIN_H 720

# define SPHERE 1
# define PLANE 2

# include "mlx.h"
# include <fcntl.h>
# include <math.h>
# include <stdio.h>
# include <time.h>
# include "../libft/libft.h"

typedef struct          s_lst
{
    double              var;
    struct s_lst        *next;
}                       t_lst;

typedef struct          s_object
{
    void                *object;
    int                 type;
    struct s_object     *next;
}                       t_object;

typedef struct          s_col
{
    double              red;
    double              green;
    double              blue;
    double              special;
}                       t_col;

typedef struct          s_vec
{
    double              x;
    double              y;
    double              z;
}                       t_vec;

typedef struct          s_light
{
    t_vec               position;
    t_col               colour;
}                       t_light;

typedef struct          s_sphere
{
    t_vec               center;
    double              radius;
    t_col               colour;
}                       t_sphere;

typedef struct          s_plane
{
    t_vec               normal;
    double              distance;
    t_col               colour;
}                       t_plane;

typedef struct          s_ray
{
    t_vec               origin;
    t_vec               direction;
}                       t_ray;

typedef struct          s_cam
{
    t_vec               cam_pos;
    t_vec               cam_dir;
    t_vec               cam_down;
    t_vec               cam_right;
}                       t_cam;

typedef struct          s_mx
{
    void                *mptr;
    void                *wptr;
    void                *iptr;
    int                 *rt;
    int                 bpp;
    int                 size;
    int                 end;
    int                 size_x;
    int                 size_y;
}                       t_mx;

int     ft_open(char *str, t_mx *v);

# endif