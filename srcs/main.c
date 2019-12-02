/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <marvin@42.fr>                    +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/11/29 00:30:45 by ahkhilad          #+#    #+#             */
/*   Updated: 2019/11/29 00:30:52 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

int     key_press(int keycode, void *p)
{
    t_mx    *r;

    r = (t_mx *)p;
    if (keycode == 53)
        exit(0);
    return (0);
}

double  magnitude(t_vec v)
{
    t_vec   c;
    c.x = v.x * v.x;
    c.y = v.y * v.y;
    c.z = v.z * v.z;
    return (sqrt(c.x + c.y + c.z));
}

t_vec  normalize(t_vec v)
{
    double  magnitude;

    magnitude = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    v.x /= magnitude;
    v.y /= magnitude;
    v.z /= magnitude;
    return (v);
}

t_vec   negative(t_vec v)
{
    v.x *= -1;
    v.y *= -1;
    v.z *= -1;
    return (v);
}

double  dotproduct(t_vec a, t_vec b)
{
    return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
}

t_vec   crossproduct(t_vec a, t_vec b)
{
    t_vec   c;

    c.x = (a.y * b.z) - (a.z * b.y);
    c.y = (a.z * b.x) - (a.x * b.z);
    c.z = (a.x * b.y) - (a.y * b.x);
    return (c);
}

t_vec   vectoradd(t_vec a, t_vec b)
{
    t_vec   c;

    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return (c);
}

t_vec   vectormulti(double scalar)
{
    t_vec   c;

    c.x *= scalar;
    c.y *= scalar;
    c.z *= scalar;
    return (c);
}

double  intersections(t_ray *r, t_plane *p)
{
    t_vec   r_dir;
    double  a;
    double  b;

    r_dir.x = r->direction.x;
    r_dir.y = r->direction.y;
    r_dir.z = r->direction.z;
    a = dotproduct(r_dir, p->normal);
    if (a == 0)
        // ray is parallel to the plane
        return (-1);
    else
    {
        p->normal = negative(vectormulti(p->distance));
        b = dotproduct(p->normal, vectoradd(r->origin, p->normal));
        return (-1 * (b / a));
    }
    return (0);
}

t_vec   normal_vec_at(t_vec point)
{
    t_plane     *p;

    p->normal.x = point.x;
    p->normal.y = point.y;
    p->normal.z = point.z;
    return (p->normal);
}

// t_vec   n_vec_at(t_vec point)
// {

// }

// ***********************************************************

// t_cam   camera(t_cam v)
// {
//     t_cam   n;

//     n.cam_pos = v.cam_pos;
//     n.cam_dir = v.cam_dir;
//     n.cam_right = v.cam_right;
//     n.cam_down = v.cam_down;
//     return (n);
// }

int     main()
{
    // t_ray   s;

    // s.r_src.x = WIN_W / 2;
    // s.r_src.y = WIN_H / 2;
    // s.r_dir.x = WIN_W / 2;
    // s.r_dir.y = (WIN_H / 4) * 3;

    // normilizing vectors;
    // s.src_length = sqrt(pow(s.r_src.x, 2) + pow(s.r_src.y, 2));
    // s.dir_length = sqrt(pow(s.r_dir.x, 2) + pow(s.r_dir.y, 2));
    // s.normalize_src = pow((s.r_src.x / s.src_length), 2) + pow((s.r_src.y / s.src_length), 2);
    // s.normalize_dir = pow((s.r_dir.x / s.dir_length), 2) + pow((s.r_dir.y / s.dir_length), 2);
    // ft_putnbr(s.normalize_src);
    // ft_putchar('\n');
    // ft_putnbr(s.normalize_dir);
    //if (s.normalize_src == 1)
        //vector is normalized
    //else
        //vector is not normalized
    //if (s.normalize_dir == 1)
        //vector is normalized
    //else
        //vector is not normalized
    t_mx        r;
    t_cam       d;
    t_cam       n;
    t_vec       o;
    t_vec       a;
    t_vec       b;
    t_vec       c;
    t_vec       look_at;
    t_vec       diff_btw;
    t_vec       light_pos;
    t_vec       cam_ray_src;
    t_vec       cam_ray_dir;
    t_col       w_light;
    t_col       green_c;
    t_col       grey_c;
    t_col       black_c;
    t_col       maroon_c;
    t_ray       cam_ray;
    t_light     light;
    t_sphere    sphere;
    t_plane     plane;
    int         i;
    int         j;
    int         color;
    double      asp_ratio;
    double      x;
    double      y;

    r.mptr = mlx_init();
    r.wptr = mlx_new_window(r.mptr, WIN_W, WIN_H, "RTV1");
    r.iptr = mlx_new_image(r.mptr, WIN_W, WIN_H);
    r.rt = (int *)mlx_get_data_addr(r.iptr, &r.bpp, &r.size, &r.end);
    i = 0;
    asp_ratio = (double)WIN_W / (double)WIN_H;
    o.x = 0;
    o.y = 0;
    o.z = 0;
    a.x = 1;
    a.y = 0;
    a.z = 0;
    b.x = 0;
    b.y = 1;
    b.z = 0;
    c.x = 0;
    c.y = 0;
    c.z = 1;
    light_pos.x = -7;
    light_pos.y = 10;
    light_pos.z = -10;

    // light------------
    
    w_light.red = 1.0;
    w_light.green = 1.0;
    w_light.blue = 1.0;
    w_light.special = 0;

    // colors-----------

    green_c.red = 0.5;
    green_c.green = 1.0;
    green_c.blue = 0.5;
    green_c.special = 0.3;
    grey_c.red = 0.5;
    grey_c.green = 0.5;
    grey_c.blue = 0.5;
    grey_c.special = 0;
    black_c.red = 0;
    black_c.green = 0;
    black_c.blue = 0;
    black_c.special = 0;
    maroon_c.red = 0.5;
    maroon_c.green = 0.25;
    maroon_c.blue = 0.25;
    maroon_c.special = 0;
    light.position.x = light_pos.x;
    light.position.y = light_pos.y;
    light.position.z = light_pos.z;
    light.colour.red = w_light.red;
    light.colour.green = w_light.green;
    light.colour.blue = w_light.blue;
    light.colour.special = w_light.special;

    // objects----------

    sphere.center.x = o.x;
    sphere.center.y = o.y;
    sphere.center.z = o.z;
    sphere.radius = 1.0;
    sphere.colour.red = green_c.red;
    sphere.colour.green = green_c.green;
    sphere.colour.blue = green_c.blue;
    sphere.colour.special = green_c.special;
    plane.normal.x = b.x;
    plane.normal.y = b.y;
    plane.normal.z = b.z;
    plane.distance = -1;
    plane.colour.red = maroon_c.red;
    plane.colour.green = maroon_c.green;
    plane.colour.blue = maroon_c.blue;
    plane.colour.special = maroon_c.special;

    // -----------------

    look_at.x = 0;
    look_at.y = 0;
    look_at.z = 0;
    diff_btw.x = d.cam_pos.x - look_at.x;
    diff_btw.y = d.cam_pos.y - look_at.y;
    diff_btw.z = d.cam_pos.z - look_at.z;
    d.cam_pos.x = 3;
    d.cam_pos.y = 1.5;
    d.cam_pos.z = -4;
    d.cam_dir.x = 0;
    d.cam_dir.y = 0;
    d.cam_dir.z = 1;
    d.cam_right.x = 0;
    d.cam_right.y = 0;
    d.cam_right.z = 0;
    d.cam_down.x = 0;
    d.cam_down.y = 0;
    d.cam_down.z = 0;
    d.cam_dir = normalize(negative(diff_btw));
    d.cam_right = normalize(crossproduct(b, d.cam_dir));
    d.cam_down = crossproduct(d.cam_right, d.cam_dir);
    n = d;
    while (i < WIN_W)
    {
        j = 0;
        while(j < WIN_H)
        {
            color = j * WIN_W + i;
            // no anti-aliasing
            if (WIN_W > WIN_H)
            {
                // the image is wider than it's taller
                x = (((i + 0.5) / WIN_W) * asp_ratio) - (((WIN_W - WIN_H) / (double)WIN_H) / 2);
                y = ((WIN_H - j) + 0.5) / WIN_H;
            }
            else if (WIN_H > WIN_W)
            {
                // the image is taller than it's wider
                x = ((i + 0.5) / WIN_W);
                y = ((((WIN_H - j) + 0.5) / WIN_H) / asp_ratio) - (((WIN_H - WIN_W) / (double)WIN_W) / 2);
            }
            else
            {
                // the image is square
                x = (i + 0.5) / WIN_W;
                y = ((WIN_H - j) + 0.5) / WIN_H;
            }
            cam_ray_src = n.cam_pos;
            n.cam_right = vectormulti(x - 0.5);
            n.cam_down = vectormulti(y - 0.5);
            cam_ray_dir = normalize(vectoradd(n.cam_dir, vectoradd(n.cam_right, n.cam_down)));
            cam_ray.origin = cam_ray_src;
            cam_ray.direction = cam_ray_dir;

            // ****************************************************

            if ((i >= (WIN_W / 4) && i <= (WIN_W - (WIN_W / 4))) && (j >= (WIN_H / 4) && j <= (WIN_H - (WIN_H / 4))))
                r.rt[color] = 0xff0055;
            else
                r.rt[color] = 0;
            j++;
        }
        i++;
    }
    mlx_put_image_to_window(r.mptr, r.wptr, r.iptr, 0, 0);
    mlx_hook(r.wptr, 2, 0, key_press, (void *)&r);
    mlx_loop(r.mptr);
    return (0);
}
