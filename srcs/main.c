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

int         key_press(int keycode, void *p)
{
    t_mx    *r;

    r = (t_mx *)p;
    if (keycode == 53)
        exit(0);
    return (0);
}

double      magnitude(t_vec v)
{
    t_vec   c;
    
    c.x = v.x * v.x;
    c.y = v.y * v.y;
    c.z = v.z * v.z;
    return (sqrt(c.x + c.y + c.z));
}

t_vec       normalize(t_vec v)
{
    double  magnitude;

    magnitude = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    v.x /= magnitude;
    v.y /= magnitude;
    v.z /= magnitude;
    return (v);
}

t_vec       negative(t_vec v)
{
    v.x *= -1;
    v.y *= -1;
    v.z *= -1;
    return (v);
}

double      dotproduct(t_vec a, t_vec b)
{
    return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
}

t_vec       crossproduct(t_vec a, t_vec b)
{
    t_vec   c;

    c.x = (a.y * b.z) - (a.z * b.y);
    c.y = (a.z * b.x) - (a.x * b.z);
    c.z = (a.x * b.y) - (a.y * b.x);
    return (c);
}

t_vec       vectoradd(t_vec a, t_vec b)
{
    t_vec   c;

    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return (c);
}

t_vec       vectormulti(double scalar)
{
    t_vec   c;

    c.x *= scalar;
    c.y *= scalar;
    c.z *= scalar;
    return (c);
}

double      find_intersection(t_ray *r, t_plane *p)
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

t_vec       ft_normal_at(t_vec point, t_sphere sphere)
{
    t_vec   normal_vect;

    normal_vect = normalize(negative(vectoradd(point, sphere.center)));
    return (normal_vect);
}

double      find_intersection_sphere(t_ray ray, t_sphere sphere)
{
    t_vec   ray_origin;
    double  ray_origin_x;
    double  ray_origin_y;
    double  ray_origin_z;
    t_vec   ray_direction;
    double  ray_direction_x;
    double  ray_direction_y;
    double  ray_direction_z;
    t_vec   sphere_center;
    double  sphere_center_x;
    double  sphere_center_y;
    double  sphere_center_z;
    double  a;
    double  b;
    double  c;
    double  discriminant;
    double  root_1;
    double  root_2;

    ray_origin = ray.origin;
    ray_origin_x = ray_origin.x;
    ray_origin_y = ray_origin.y;
    ray_origin_z = ray_origin.z;

    ray_direction = ray.direction;
    ray_direction_x = ray_direction.x;
    ray_direction_y = ray_direction.y;
    ray_direction_z = ray_direction.z;

    sphere_center = sphere.center;
    sphere_center_x = sphere_center.x;
    sphere_center_y = sphere_center.y;
    sphere_center_z = sphere_center.z;
    a = 1;  // normalized
    b = (2 * (ray_origin_x - sphere_center_x) * ray_direction_x)\
    + (2 * (ray_origin_y - sphere_center_y) * ray_direction_y)\
    + (2 * (ray_origin_z - sphere_center_z) * ray_direction_z);
    c = pow((ray_origin_x - sphere_center_x), 2)\
    + pow((ray_origin_y - sphere_center_y), 2)\
    + pow((ray_origin_z - sphere_center_z), 2) - (sphere.radius * sphere.radius);
    discriminant = (b * b) - (4 * c);
    if (discriminant > 0)
    {
        // the ray intersects the sphere
        // first root
        root_1 = ((-1 * b - sqrt(discriminant)) / 2) - 0.000001;
        if (root_1 > 0)
        {
            // the first root is the smallest positive root
            return (root_1);
        }
        else
        {
            // the second root is the smallest positive root
            root_2 = ((sqrt(discriminant) - b) / 2) - 0.000001;
            return (root_2);
        }    
    }
    else
    {
        // the ray missed the sphere
        return (-1);
    }
}

t_vec       normal_vec_at(t_vec point)
{
    t_plane     *p;

    p->normal.x = point.x;
    p->normal.y = point.y;
    p->normal.z = point.z;
    return (p->normal);
}

void        ft_add_node_at_last(t_lst **head, t_lst *new)
{
	t_lst *tmp;

	tmp = *head;
    if (!tmp)
        *head = new;
    else
		ft_add_node_at_last(&(tmp->next), new);
}

void        ft_add_node_at_last_obj(t_object **head, t_object *new)
{
	t_object *tmp;

	tmp = *head;
    if (!tmp)
        *head = new;
    else
		ft_add_node_at_last_obj(&(tmp->next), new);
}

t_object    *ft_malloc_node_obj(void *obj, int type)
{
    t_object    *node;

    node = (t_object*)malloc(sizeof(t_object));
    node->object = obj;
    node->type = type;
    node->next = NULL;
    return(node);
}

t_lst    *ft_malloc_node_lst(double var)
{
    t_lst    *node;

    node = (t_lst*)malloc(sizeof(t_lst));
    node->var = var;
    node->next = NULL;
    return(node);
}



// t_sphere *ft_get_data_sphere()
// {
//     char    *line;
//     t_sphere *sp;
//     sp = malloc(sizeof(t_sphere));
//     while (get_next_line(line))
//     {
//         sp->radius = 5;
//         sp->center = (t_vec){10,0,1};
//         free(line);
//     }
// }

// void ft_parse_file(t_list *obj, char *line)
// {
//     void *content;
//     if (strcmp(line, 'sphere'))
//     {
//         content = ft_get_data_sphere();
//         obj = ft_lstnew(content, SPHERE);
//     }
//     // else return error 
// }

// void    make_list()
// {
//     t_list *head = NULL;
//     t_list *obj;
//     char    *line;
//     while (get_next_line(line))
//     {
//         ft_parse_file(&obj, line);
//             // return ('error')
//         ft_lstadd(&head, obj);
//         obj = obj->next;
//         free(line);
//     }
//     exit(1);
// }


int         main(void)
{
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
    t_object    *head;
    int         index;
    int         object_index;
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

    // light_color------------
    
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

    // light------------

    light.position = light_pos;
    light.colour = w_light;

    // objects----------

    sphere.center = o;
    sphere.radius = 1.0;
    sphere.colour = green_c;
    plane.normal = b;
    plane.distance = -1;
    plane.colour = maroon_c;

    // -----------------

    ft_add_node_at_last_obj(&head, ft_malloc_node_obj(&sphere, SPHERE));
    ft_add_node_at_last_obj(&head, ft_malloc_node_obj(&plane, PLANE));

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
            
            index = 0;

            // while (index < 2)
            // {
            //     // ************************** something here
            //     [index] = ;
            //     index += 1;
            // }
            // object_index = winning_object_index();

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
