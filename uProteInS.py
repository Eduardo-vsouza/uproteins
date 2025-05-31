#!/usr/bin/python3

# Copyright © 2021-2025 Eduardo Vieira de Souza
# Copyright © 2021-2025 Adriana Canedo
# Copyright © 2021-2025 Cristiano Valim Bizarro
# Copyright © 2025 Bruno Maestri A Becker
#
# This file is part of uProteInS.
#
# uProteInS is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# uProteInS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# uProteInS. If not, see <https://www.gnu.org/licenses/>.


import typing as t

from src import cli, main


def uproteins(args: t.Optional[t.Sequence[str]] = None):
    try:
        from pyfiglet import Figlet
        f = Figlet()
        print(f.renderText('uProteInS'))
    except ModuleNotFoundError:
        print('uProteInS')

    parser, subparsers = cli.get_parsers()
    args_namespace = parser.parse_args(args)
    main.run_workflow(args_namespace, subparsers[args_namespace.mode])


if __name__ == "__main__":
    uproteins()
