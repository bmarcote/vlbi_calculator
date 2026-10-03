"""Shared argparse helpers for renamed CLI options (stdlib only, so `planobs -h` stays fast).

Two kinds of legacy spellings are supported, uniformly across every `planobs` subcommand:
- Deprecated aliases (`add_deprecated_alias`): hidden from the help, still work, store to the new option's
  dest and print a one-line deprecation warning to stderr.
- Removed options (`add_removed_option`): hidden from the help; using them aborts the parsing with an
  explicit error naming the new option. Used when the old spelling now means something else elsewhere.
"""
import sys
import argparse


def deprecation_message(old_flags: list[str] | tuple[str, ...], new_option: str) -> str:
    """Return the one-line deprecation warning, e.g. "Warning: '-t1/--starttime' is deprecated; use '-e/--epoch'."."""
    return f"Warning: '{'/'.join(old_flags)}' is deprecated; use '{new_option}'."


class DeprecatedOptionAction(argparse.Action):
    """Hidden alias of a renamed option: stores the value into the new option's dest and warns on stderr.

    Extra constructor keyword: new_option (str), the new spelling shown in the warning (e.g. '-e/--epoch').
    The default is SUPPRESS so the alias never overwrites the default set by the new option.
    """

    def __init__(self, option_strings, dest, new_option: str, **kwargs):
        kwargs['help'] = argparse.SUPPRESS
        kwargs['default'] = argparse.SUPPRESS
        super().__init__(option_strings, dest, **kwargs)
        self.new_option = new_option

    def __call__(self, parser, namespace, values, option_string=None):
        print(deprecation_message(self.option_strings, self.new_option), file=sys.stderr)
        setattr(namespace, self.dest, values)


class RemovedOptionAction(argparse.Action):
    """Hidden former option whose use is an error (its spelling now has another meaning in other commands).

    Extra constructor keywords: new_option (str) and reason (str), both shown in the error message.
    Accepts an optional value (nargs='?') so both '-t' and '-t VALUE' reach the explicit error.
    """

    def __init__(self, option_strings, dest, new_option: str, reason: str, **kwargs):
        kwargs['help'] = argparse.SUPPRESS
        kwargs['default'] = argparse.SUPPRESS
        kwargs['nargs'] = '?'
        super().__init__(option_strings, dest, **kwargs)
        self.new_option = new_option
        self.reason = reason

    def __call__(self, parser, namespace, values, option_string=None):
        parser.error(f"'{option_string}' is no longer accepted here ({self.reason}); use '{self.new_option}'.")


def add_deprecated_alias(parser, *flags: str, dest: str, new_option: str, **kwargs) -> None:
    """Add hidden deprecated spellings `flags` that store into `dest` (the new option's dest) and warn.

    Parameters
    ----------
    parser : argparse.ArgumentParser or argparse._ArgumentGroup
        Where to add the alias.
    *flags : str
        Old option strings, e.g. '-t1', '--starttime'.
    dest : str
        Destination attribute of the new option (e.g. 'epoch').
    new_option : str
        New spelling shown in the warning (e.g. '-e/--epoch').
    **kwargs
        Extra add_argument keywords that must match the new option (type, nargs).
    """
    parser.add_argument(*flags, dest=dest, action=DeprecatedOptionAction, new_option=new_option, **kwargs)


def add_removed_option(parser, flag: str, new_option: str, reason: str) -> None:
    """Add a hidden former option `flag` whose use stops parsing with an error pointing to `new_option`.

    Parameters
    ----------
    parser : argparse.ArgumentParser or argparse._ArgumentGroup
        Where to add the option.
    flag : str
        Former option string (e.g. '-t').
    new_option : str
        New spelling the user must use instead (e.g. '-e/--epoch').
    reason : str
        Short explanation shown in the error (e.g. "'-t' means '--target' in other commands").
    """
    parser.add_argument(flag, dest=f"_removed_{flag.lstrip('-')}", action=RemovedOptionAction,
                        new_option=new_option, reason=reason)
