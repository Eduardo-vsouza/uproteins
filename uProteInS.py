#!/usr/bin/python3

from src import cli, main


def uproteins():
    try:
        from pyfiglet import Figlet
        f = Figlet()
        print(f.renderText('uProteInS'))
    except ModuleNotFoundError:
        pass

    parser, subparsers = cli.get_parsers()
    args = parser.parse_args()
    main.run_workflow(args, subparsers[args.mode])


if __name__ == "__main__":
    uproteins()
