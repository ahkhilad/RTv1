/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_atof.c                                          :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahkhilad <ahkhilad@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/03/27 12:03:27 by ahkhilad          #+#    #+#             */
/*   Updated: 2020/06/21 21:12:13 by ahkhilad         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

double		ft_atof(const char *str)
{
	double			whole_part;
	double			decimal_part;
	char			*point;
	int				i;

	if (!(point = ft_strchr(str, '.')) || (point[1] && !ft_isdigit(point[1])))
		return (ft_atoi(str));
	whole_part = ft_atoi(str);
	decimal_part = ft_atoi(&point[1]);
	i = 1;
	while (point[i])
	{
		if (!ft_isdigit(point[i]))
			break ;
		i++;
	}
	while (i > 1)
	{
		decimal_part = decimal_part / 10;
		i--;
	}
	return (whole_part + ((whole_part >= 0) ? decimal_part : -decimal_part));
}
