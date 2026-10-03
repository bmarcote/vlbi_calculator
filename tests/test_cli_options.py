"""Tests for the unified CLI option names across all planobs subcommands and their legacy spellings."""
import argparse
import sys

import pytest

from vlbiplanobs import cli

EPOCH = '2026-11-10 10:00'
DEPRECATION_EPOCH = "'-t1/--starttime' is deprecated; use '-e/--epoch'"


def _parser(add_arguments) -> argparse.ArgumentParser:
    """Return a bare parser populated by the given cli add_*_arguments function."""
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    return parser


def _run_cli(monkeypatch, argv: list[str], handler_name: str) -> argparse.Namespace:
    """Run cli.cli() with `argv`, capturing the Namespace passed to the (patched) handler `handler_name`."""
    captured = {}
    monkeypatch.setattr(cli, handler_name, lambda args: captured.setdefault('args', args))
    monkeypatch.setattr(sys, 'argv', ['planobs'] + argv)
    cli.cli()
    return captured['args']


# New option names: each subcommand parses to the expected dest.

def test_observe_new_flags_parse_to_expected_dests():
    args = _parser(cli.add_observation_arguments).parse_args(
        ['-t', 'J1230+1223', 'M87', '-b', '18cm', '-n', 'EVN', '-e', EPOCH, '-d', '8'])
    assert args.targets == ['J1230+1223', 'M87'] and args.epoch == EPOCH and args.duration == 8.0
    assert args.network == ['EVN'] and args.band == '18cm'
    args = _parser(cli.add_observation_arguments).parse_args(['--target', 'M87', '--epoch', EPOCH, '--network', 'EVN'])
    assert args.targets == ['M87'] and args.epoch == EPOCH and args.network == ['EVN']


def test_fringefinders_new_flags_parse_to_expected_dests():
    args = _parser(cli.add_fringe_finder_arguments).parse_args(['-n', 'EVN', '-e', EPOCH, '-d', '8', '-l', '5'])
    assert args.network == ['EVN'] and args.epoch == EPOCH and args.duration == 8.0 and args.max_lines == 5
    args = _parser(cli.add_fringe_finder_arguments).parse_args(['--epoch', EPOCH, '-d', '1', '--max-lines', '3'])
    assert args.epoch == EPOCH and args.max_lines == 3


def test_phasecals_new_flags_parse_to_expected_dests():
    args = _parser(cli.add_phase_cal_arguments).parse_args(['J1230+1223', '-b', '6cm', '-l', '5'])
    assert args.target == 'J1230+1223' and args.max_lines == 5 and args.band == '6cm'
    args = _parser(cli.add_phase_cal_arguments).parse_args(['-t', 'J1230+1223', '--max-lines', '2'])
    assert args.target is None and args.target_option == 'J1230+1223' and args.max_lines == 2
    assert _parser(cli.add_phase_cal_arguments).parse_args(['M87']).max_lines is None


def test_legacy_default_mode_uses_new_flags(monkeypatch):
    args = _run_cli(monkeypatch, ['-t', 'M87', '-b', '18cm', '-n', 'EVN', '-e', EPOCH, '-d', '8'],
                    'handle_observation_command')
    assert args.command == 'observe' and args.targets == ['M87'] and args.epoch == EPOCH and args.duration == 8.0


def test_subcommands_use_new_flags(monkeypatch):
    args = _run_cli(monkeypatch, ['observe', '-t', 'M87', '-b', '18cm', '-n', 'EVN', '-e', EPOCH, '-d', '8'],
                    'handle_observation_command')
    assert args.targets == ['M87'] and args.epoch == EPOCH
    args = _run_cli(monkeypatch, ['fringefinders', '-n', 'EVN', '-e', EPOCH, '-d', '8', '-l', '5'],
                    'handle_fringe_finder_command')
    assert args.epoch == EPOCH and args.max_lines == 5
    args = _run_cli(monkeypatch, ['phasecals', '-t', 'M87', '-l', '4'], 'handle_phase_cal_command')
    assert args.target_option == 'M87' and args.max_lines == 4


# Deprecated aliases: still work, store into the new dest and warn on stderr.

@pytest.mark.parametrize('flag', ['-t1', '--starttime'])
def test_observe_starttime_alias_warns(flag, capsys):
    args = _parser(cli.add_observation_arguments).parse_args([flag, EPOCH, '-b', '18cm'])
    assert args.epoch == EPOCH and not hasattr(args, 'starttime')
    assert DEPRECATION_EPOCH in capsys.readouterr().err


def test_observe_targets_alias_warns(capsys):
    args = _parser(cli.add_observation_arguments).parse_args(['--targets', 'M87', 'J1230+1223'])
    assert args.targets == ['M87', 'J1230+1223']
    assert "'--targets' is deprecated; use '-t/--target'" in capsys.readouterr().err


def test_new_flags_do_not_warn(capsys):
    _parser(cli.add_observation_arguments).parse_args(['-t', 'M87', '-e', EPOCH])
    assert capsys.readouterr().err == ''


@pytest.mark.parametrize('flag', ['-t1', '--starttime'])
def test_fringefinders_starttime_alias_warns(flag, capsys):
    args = _parser(cli.add_fringe_finder_arguments).parse_args(['-n', 'EVN', flag, EPOCH, '-d', '4'])
    assert args.epoch == EPOCH
    assert DEPRECATION_EPOCH in capsys.readouterr().err


def test_phasecals_n_sources_alias_warns(capsys):
    args = _parser(cli.add_phase_cal_arguments).parse_args(['-t', 'M87', '--n-sources', '3'])
    assert args.max_lines == 3
    assert "'--n-sources' is deprecated; use '-l/--max-lines'" in capsys.readouterr().err


def test_legacy_default_mode_accepts_t1_alias(monkeypatch, capsys):
    args = _run_cli(monkeypatch, ['-t', 'M87', '-b', '18cm', '-n', 'EVN', '-t1', EPOCH, '-d', '8'],
                    'handle_observation_command')
    assert args.epoch == EPOCH and args.targets == ['M87']
    assert DEPRECATION_EPOCH in capsys.readouterr().err


# Removed spellings whose letter now means something else: explicit error.

@pytest.mark.parametrize('argv', [['-t', EPOCH], ['-t']])
def test_fringefinders_t_is_an_error(argv, capsys):
    with pytest.raises(SystemExit) as exc:
        _parser(cli.add_fringe_finder_arguments).parse_args(['-n', 'EVN', '-d', '4'] + argv)
    assert exc.value.code == 2
    assert "'-t' is no longer accepted" in capsys.readouterr().err


def test_phasecals_n_is_an_error(capsys):
    with pytest.raises(SystemExit) as exc:
        _parser(cli.add_phase_cal_arguments).parse_args(['M87', '-n', '5'])
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "'-n' is no longer accepted" in err and "-l/--max-lines" in err


def test_fringefinders_missing_epoch_is_an_error(capsys):
    args = _parser(cli.add_fringe_finder_arguments).parse_args(['-n', 'EVN', '-d', '4'])
    with pytest.raises(SystemExit) as exc:
        cli.handle_fringe_finder_command(args)
    assert exc.value.code == 1
    assert '-e/--epoch' in capsys.readouterr().out


def test_deprecated_options_hidden_from_help():
    for add_arguments in (cli.add_observation_arguments, cli.add_fringe_finder_arguments,
                          cli.add_phase_cal_arguments):
        help_text = _parser(add_arguments).format_help()
        for old in ('-t1', '--starttime', '--targets', '--n-sources'):
            assert old not in help_text


def test_phasecals_rfc_catalog_and_deprecated_catalog_file(capsys):
    """`--rfc-catalog` is the option name; `--catalog-file` still works but warns."""
    import argparse
    from vlbiplanobs import cli
    parser = argparse.ArgumentParser()
    cli.add_phase_cal_arguments(parser)
    assert parser.parse_args(['J1230+1223', '--rfc-catalog', 'x.txt']).rfc_catalog == 'x.txt'
    assert parser.parse_args(['J1230+1223', '--catalog-file', 'y.txt']).rfc_catalog == 'y.txt'
    assert '--catalog-file' in capsys.readouterr().err
