#!/bin/bash

scene_element_options() {

    tput setaf 9;
    echo -e "allowed elements :
    => sphere
    => plane
    => light
    => ray
    => camera"

    read option

    if [[ ( $option == "sphere" || $option == "plane" || $option == "light" || $option == "ray" || $option == "camera" ) ]]
    then
    echo "$option" >> ./"$file_name_str".conf
        if [ $option == "sphere" ];
        then
            echo -e "please add its center vector !"
            read -p "x = " x
            while [ true ]
            do
                if ! [[ "$x" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "x = " x
                else
                break
                fi
            done
            read -p "y = " y
            while [ true ]
            do
                if ! [[ "$y" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "y = " y
                else
                break
                fi
            done
            read -p "z = " z
            while [ true ]
            do
                if ! [[ "$z" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "z = " z
                else
                break
                fi
            done
            echo -e "\tcenter $x $y $z" >> ./"$file_name_str".conf
            echo -e "please add its radius !"
            read -p "radius = " radius
            while [ true ]
            do
                if ! [[ "$radius" =~ ^[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry positive integers only !"
                    tput setaf 9;
                    read -p "radius = " radius
                else
                break
                fi
            done
            echo -e "\tradius $radius" >> ./"$file_name_str".conf
            echo -e "please add its color in Hexadecimal !"
            read -p "color = " color
            echo -e "\tcolor $color" >> ./"$file_name_str".conf
            echo -e "." >> ./"$file_name_str".conf
        add_new_element_function
        elif [ $option == "plane" ];
        then
            echo -e "please add its normal vector"
            read -p "x = " x
            while [ true ]
            do
                if ! [[ "$x" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "x = " x
                else
                break
                fi
            done
            read -p "y = " y
            while [ true ]
            do
                if ! [[ "$y" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "y = " y
                else
                break
                fi
            done
            read -p "z = " z
            while [ true ]
            do
                if ! [[ "$z" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "z = " z
                else
                break
                fi
            done
            echo -e "\tnormal $x $y $z" >> ./"$file_name_str".conf
            echo -e "please add its distance !"
            read -p "distance = " distance
            while [ true ]
            do
                if ! [[ "$distance" =~ ^[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry positive integers only !"
                    tput setaf 9;
                    read -p "distance = " distance
                else
                break
                fi
            done
            echo -e "\tdistance $distance" >> ./"$file_name_str".conf
            echo -e "please add its color in Hexadecimal !"
            read -p "color = " color
            echo -e "\tcolor $color" >> ./"$file_name_str".conf
            echo -e "." >> ./"$file_name_str".conf
        add_new_element_function
        elif [ $option == "light" ];
        then
            echo -e "please add its position vector"
            read -p "x = " x
            while [ true ]
            do
                if ! [[ "$x" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "x = " x
                else
                break
                fi
            done
            read -p "y = " y
            while [ true ]
            do
                if ! [[ "$y" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "y = " y
                else
                break
                fi
            done
            read -p "z = " z
            while [ true ]
            do
                if ! [[ "$z" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "z = " z
                else
                break
                fi
            done
            echo -e "\tposition $x $y $z" >> ./"$file_name_str".conf
            echo -e "please add its color in Hexadecimal !"
            read -p "color = " color
            echo -e "\tcolor $color" >> ./"$file_name_str".conf
            echo -e "." >> ./"$file_name_str".conf
        add_new_element_function
        elif [ $option == "ray" ];
        then
            echo -e "please add its origin vector"
            read -p "x = " x
            while [ true ]
            do
                if ! [[ "$x" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "x = " x
                else
                break
                fi
            done
            read -p "y = " y
            while [ true ]
            do
                if ! [[ "$y" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "y = " y
                else
                break
                fi
            done
            read -p "z = " z
            while [ true ]
            do
                if ! [[ "$z" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "z = " z
                else
                break
                fi
            done
            echo -e "\torigin $x $y $z" >> ./"$file_name_str".conf
            echo -e "please add its direction vector"
            read -p "x = " x
            while [ true ]
            do
                if ! [[ "$x" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "x = " x
                else
                break
                fi
            done
            read -p "y = " y
            while [ true ]
            do
                if ! [[ "$y" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "y = " y
                else
                break
                fi
            done
            read -p "z = " z
            while [ true ]
            do
                if ! [[ "$z" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "z = " z
                else
                break
                fi
            done
            echo -e "\tdirection $x $y $z" >> ./"$file_name_str".conf
            echo -e "." >> ./"$file_name_str".conf
        add_new_element_function
        elif [ $option == "camera" ];
        then
            echo -e "please add its position vector !"
            read -p "x = " x
            while [ true ]
            do
                if ! [[ "$x" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "x = " x
                else
                break
                fi
            done
            read -p "y = " y
            while [ true ]
            do
                if ! [[ "$y" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "y = " y
                else
                break
                fi
            done
            read -p "z = " z
            while [ true ]
            do
                if ! [[ "$z" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "z = " z
                else
                break
                fi
            done
            echo -e "\tposition $x $y $z" >> ./"$file_name_str".conf
            echo -e "please add its direction vector !"
            read -p "x = " x
            while [ true ]
            do
                if ! [[ "$x" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "x = " x
                else
                break
                fi
            done
            read -p "y = " y
            while [ true ]
            do
                if ! [[ "$y" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "y = " y
                else
                break
                fi
            done
            read -p "z = " z
            while [ true ]
            do
                if ! [[ "$z" =~ ^-?[0-9]+$ ]]
                then
                    tput setaf 1;
                    echo "Sorry integers only !"
                    tput setaf 9;
                    read -p "z = " z
                else
                break
                fi
            done
            echo -e "\tdirection $x $y $z" >> ./"$file_name_str".conf
            echo -e "." >> ./"$file_name_str".conf
        add_new_element_function
        fi
    else
        tput setaf 1;
        echo "please add one of the allowed elements"
        scene_element_options
    fi
}

file_name_function() {
tput setaf 4;
    read -p "please enter the file name, otherwise press Q to exit : " file_name_str
        if [[ ( "$file_name_str" == "Q" || "$file_name_str" == "q" ) ]]
        then
            exit 0
        elif [[ ! -z "$file_name_str" ]]
        then
        tput setaf 2;
            echo "file name has been accepetd !"
        else
        tput setaf 1;
            echo "file name cannot be empty !"
            file_name_function
        fi
}

add_new_element_function() {

    tput setaf 4;
    read -p "if u want to add another element, please press C, otherwise press Q to exit : " choice
    if [[ ! -z "$choice" ]]
    then
        if [[ ( $choice == C || $choice == c ) ]];
        then
            scene_element_options
        elif [[ ( $choice == Q || $choice == q ) ]];
        then
            length=$(wc -c < "$file_name_str".conf)
            if [ "$length" -ne 0 ] && [ -z "$(tail -c -1 < "$file_name_str".conf)" ]; then
            # The file ends with a newline or null
            dd if=/dev/null of="$file_name_str".conf obs="$((length-1))" seek=1
            fi 
            exit 0
        fi
    else
        tput setaf 1;
        echo "wrong argument !"
        add_new_element_function
    fi
}

file_name_function
scene_element_options