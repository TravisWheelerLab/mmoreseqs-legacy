#!/usr/bin/env python
###############################################################################
#  - FILE: 	application.py
#  - DESC:  Interface for using MMORESEQS search.
###############################################################################

# from email.policy import default
import os
import sys
from sre_constants import FAILURE
from tokenize import Double
import click
import json
import typing
import functools

# Include mmoreseqs module path
py_mmoreseqs_module_path = os.getcwd() + "/../build"
sys.path.append(py_mmoreseqs_module_path)

# Try to import mmoreseqs
try:
    import py_mmoreseqs
    try:
        ans = py_mmoreseqs.add(3, 4)
    except Exception as err:
        print(err)
except Exception as err:
    print(err)

debug = False


def dbg_print(*args, **kwargs):
    """
    Print if debugging. 
    """
    if debug:
        print(*args, **kwargs)


def quote(my_str: str):
    """
    Enquote string.
    """
    my_str = f'"{my_str}"'
    return my_str


def type_list(my_list: list, my_type: bool):
    """ 
    Unquote list.
    """
    my_str = "("
    for i in range(len(my_list)):
        if my_type == str:
            my_str += quote(str(my_list[i]))
        else:
            my_str += str(my_list[i])
        if i < len(my_list) - 1:
            my_str += ", "
    my_str += ")"
    return my_str


def get_type(type_name: str):
    """
    Translates type string to type.
    """
    if type_name == "bool":
        return bool
    elif type_name == "str":
        return str
    elif type_name == "int":
        return int
    elif type_name == "double":
        return Double
    else:
        raise TypeError(
            "TypeName does not correspond to valid type:", type_name)


def get_types(type_names: list):
    """
    Translates list of type strings to type.
    """
    types = []
    for type_name in type_names:
        types.append(get_type(type_name))
    return types


def get_default(type_name: str):
    """
    Gets a default for type.
    """
    if type_name == "bool":
        return False
    elif type_name == "str":
        return ""
    elif type_name == "int":
        return 0
    elif type_name == "double":
        return 0.0
    else:
        raise TypeError(
            "TypeName does not correspond to valid type:", type_name)


def get_defaults(type_names: list):
    """
    Translates list of type strings to type.
    """
    types = []
    for type_name in type_names:
        types.append(get_default(type_name))
    return types


def copy_main_dict_into_temp_dict(temp_dict: dict, main_dict: dict, dict_name: str):
    for key in temp_dict:
        if key in main_dict:
            temp_dict[key] = main_dict[key]
    # cast type name string into type
    if "type" in main_dict:
        temp_dict["type_name"] = main_dict["type"]
        temp_dict["nargs"] = len(main_dict["type"])
        if type(main_dict["type"]) == list:
            temp_dict["type"] = get_types(main_dict["type"])
        else:
            temp_dict["type"] = get_type(main_dict["type"])
    # cast default to its proper type
    if "default" in main_dict:
        for i in range(len(main_dict["default"])):
            try:
                temp_dict["default"][i] = temp_dict["type"][i](
                    temp_dict["default"][i])
            except Exception as err:
                print(f"OptionDefaultCastError: {dict_name}: {err}")
                temp_dict["default"][i] = None
                print("main_dict:\n", main_dict)
                print("temp_dict:\n", temp_dict)
    return temp_dict


def add_group(group_name: str):
    @click.group(group_name)
    def cli():
        """
        MMOREseqs command suite.
        """
    return cli


def write_command(cmd_name: str, arg_names: list, opt_names: list, func, cli_group: click.Group, cli_dict: dict):
    """
    Build command for cli_group with arguments and options.
    """
    func_str = ""
    func_str += add_command(cmd_name, func, cli_group, cli_dict, True) + "\n"
    for arg_name in arg_names:
        func_str += add_argument(arg_name, cmd_name,
                                 func, cli_group, cli_dict, True) + "\n"
    for opt_name in opt_names:
        func_str += add_option(opt_name, cmd_name, func,
                               cli_group, cli_dict, True) + "\n"

    return func_str


def build_command(cmd_name: str, arg_names: list, opt_names: list, func, cli_group: click.Group, cli_dict: dict):
    """
    Build command for cli_group with arguments and options.
    """
    for opt_name in opt_names:
        func = add_option(opt_name, cmd_name, func, cli_group, cli_dict)
    for arg_name in arg_names:
        func = add_argument(arg_name, cmd_name, func, cli_group, cli_dict)
    func = add_command(cmd_name, func, cli_group, cli_dict)
    return func


def add_command(cmd_name: str, func, cli_group: click.Group,  cli_dict: dict, to_string: bool = False):
    """
    Adds command with given name as decorator to function.
    """
    cmd_dict = cli_dict["commands"][cmd_name]
    my_cmd_dict = {
        "help": None
    }
    copy_main_dict_into_temp_dict(my_cmd_dict, cmd_dict, cmd_name)

    if to_string:
        func_str = f'@{cli_group.name}.command(name={quote(cmd_name)}, help={quote(my_cmd_dict["help"])})'
        return func_str

    @cli_group.command(name=cmd_name, help=my_cmd_dict["help"])
    def wrapped(*args, **kwargs):
        dbg_print(f"wrapped command: {cmd_name}")
        return func(*args, **kwargs)
    return wrapped


def add_argument(arg_name: str, cmd_name: str, func, cli_group: click.Group, cli_dict: dict, to_string: bool = False):
    """
    Adds argument with given name as decorator to command function.
    """
    arg_dict = cli_dict["commands"][cmd_name]["arguments"][arg_name]
    my_arg_dict = {
        "type": "str",
        "help": ""
    }
    copy_main_dict_into_temp_dict(my_arg_dict, arg_dict, arg_name)

    if to_string:
        func_str = f'@click.argument({quote(arg_name)}, type={my_arg_dict["type_name"]})'
        return func_str

    @click.argument(arg_name, type=my_arg_dict["type"])
    def wrapped(*args, **kwargs):
        dbg_print(f"wrapped argument: {arg_name}")
        return func(*args, **kwargs)
    return wrapped


def add_option(opt_name: str, cmd_name: str, func, cli_group: click.Group, cli_dict: dict, to_string: bool = False):
    """
    Adds option with given name as decorator to command function.
    """
    opt_dict = cli_dict["options"][opt_name]
    my_opt_dict = {
        "type": ["str"],
        "default": [""],
        "help": "",
        "hidden": False,
    }
    copy_main_dict_into_temp_dict(my_opt_dict, opt_dict, opt_name)

    if to_string:
        func_str = f'@click.option({quote(opt_name)}, type={type_list(my_opt_dict["type_name"], bool)}, default={type_list(my_opt_dict["default"], my_opt_dict["type"][0])}, help={quote(my_opt_dict["help"])}, hidden={my_opt_dict["hidden"]}, nargs={my_opt_dict["nargs"]})'
        return func_str

    @click.option(opt_name, type=tuple(my_opt_dict["type"]), default=tuple(my_opt_dict["default"]), help=my_opt_dict["help"], hidden=my_opt_dict["hidden"], nargs=my_opt_dict["nargs"])
    def wrapped(*args, **kwargs):
        dbg_print(f"wrapped option: {opt_name}")
        return func(*args, **kwargs)
    return wrapped


def build_click_interface_from_json(json_filepath):
    """
    Build a click interface from JSON input file.
    """
    cli = add_group("cli")

    click_json = json.load(open(json_filepath))
    try:
        command_dict = click_json["commands"]
        option_dict = click_json["options"]
    except:
        print("KeyError: bad json file -- missing parent 'command' and 'option' keys.")
    for cmd_name in command_dict.keys():
        arguments = []
        options = []
        for arg_name in command_dict[cmd_name]["arguments"]:
            arguments.append(arg_name)
        for opt_name in option_dict.keys():
            try:
                if (opt_name.startswith("//")):
                    continue
                contains_option = False
                if "all" in option_dict[opt_name]["commands"]:
                    contains_option = True
                elif cmd_name in option_dict[opt_name]["commands"]:
                    contains_option = True
                if contains_option:
                    options.append(opt_name)
            except:
                print(
                    f"KeyError: bad json file -- option '{opt_name}' missing 'commands' field.")

        def base_func():
            dbg_print("This is the base function: ", cmd_name)
            pass
        build_command(cmd_name, arguments, options,
                      base_func, cli, click_json)
        command_str = write_command(cmd_name, arguments,
                                    options, base_func, cli, click_json)
        print(f"COMMAND:\n{command_str}")
    return cli


if __name__ == '__main__':
    json_filepath = os.getcwd() + "/application_options.json"
    cli = build_click_interface_from_json(json_filepath)
    cli()
