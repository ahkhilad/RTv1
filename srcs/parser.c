/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   parser.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <marvin@42.fr>                    +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/01/04 21:42:36 by ahkhilad          #+#    #+#             */
/*   Updated: 2020/01/04 21:42:44 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

int		ft_open(char *str, t_mx *v)
{
	int		fd;
	int		i;
    int     j;
	int		k;
	int		rd;
	int		x_size;
	int		y_size;
    char    *sphere;
    char    *plane;
	char	*tmp;
	char	*line;
	char	**blocks;
    char    ***obj_spec;
	char	buff[BUFF_SIZE + 1];

	fd = open(str, O_RDONLY);
	if (read(fd, buff, 0) < 0)
		return (0);
	line = ft_strnew(0);
	while ((rd = read(fd, buff, BUFF_SIZE)) > 0)
	{
		buff[rd] = '\0';
		tmp = line;
		line = ft_strjoin(line, buff);
		free(tmp);
	}
    //ft_putstr(line);
	if (ft_strlen(line) == 0)
		return 0;
	tmp = line;
	blocks = ft_strsplit(line, '.');
	free(tmp);
    //ft_putstr(blocks[1]);
	x_size = 0;
	while (blocks[x_size] != NULL)
		x_size++;
    i = 0;
    obj_spec = (char ***)malloc(sizeof(char**) * x_size);
    while (i < x_size && blocks[i])
    {
        k = 0;
        while (blocks[i][k] && blocks[i][k] != '.')
        {
            obj_spec[i] = ft_strsplit(blocks[i], '\n');
            k++;
        }
        j = k;
        i++;
    }
    sphere = "sphere";
    plane = "plane";
    int l = 0;
    int n;
    while (l < x_size && obj_spec[l])
    {
        n = 0;
        while (n < j && obj_spec[l][n])
        {
            ft_putendl(obj_spec[l][n]);   
            // if (ft_strcmp(obj_spec[l][n], sphere) == 0)
            // {
            //     ft_putendl("this is the sphere");
            //     break;
            // }
            // if (ft_strcmp(obj_spec[l][n],plane) == 0)
            // {
            //     ft_putendl("this is the plane");
            //     break;
            // }
            // else
            //     return (0);
            n++;
        }
        l++;
    }
	return (1);
}
