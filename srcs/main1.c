/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main1.c                                            :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <marvin@42.fr>                    +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/01/04 21:42:11 by ahkhilad          #+#    #+#             */
/*   Updated: 2020/01/04 21:42:23 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

void	ft_usage(void)
{
	ft_putstr("Usage: binary file [a valid map].\n");
	exit(0);
}

void	ft_destroy(t_mx *v)
{
	//mlx_destroy_image(v->mptr, v->iptr);
	mlx_clear_window(v->mptr, v->wptr);
	mlx_destroy_window(v->mptr, v->wptr);
}

int		key_press(int keycode, void *p)
{
	t_mx	*v;

	v = (t_mx *)p;
	if (keycode == 53)
	{
		ft_destroy(v);
		exit(0);
	}
	return(0);
}

int		red_button(void *p)
{
	t_mx	*v;

	v = (t_mx *)p;
	ft_destroy(v);
	exit(0);
}

int     main(int ac, char **av)
{
    t_mx    v;

    if (ac == 2)
    {
        if (av[1])
        {
            if (!ft_open(av[1], &v))
            {
                ft_putstr("error! please try a valid map.\n");
				exit(0);
            }
            else
			{
				v.mptr = mlx_init();
				v.wptr = mlx_new_window(v.mptr, WIN_W, WIN_H, "RTv1");
				mlx_hook(v.wptr, 2, 0, key_press, (void *)&v);
				mlx_hook(v.wptr, 17, 0, red_button, (void *)&v);
				mlx_loop(v.mptr);
			}
        }
    }
    else
        ft_usage();
    return (0);
}
