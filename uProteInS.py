#!/usr/bin/python3
import typing as t

from src import cli, main


def uproteins(args: t.Optional[t.Sequence[str]] = None):
    try:
        from pyfiglet import Figlet
        f = Figlet()
        print(f.renderText('uProteInS'))
    except ModuleNotFoundError:
        pass

    parser, subparsers = cli.get_parsers()
    args = parser.parse_args(args)
    main.run_workflow(args, subparsers[args.mode])


if __name__ == "__main__":
    uproteins()
