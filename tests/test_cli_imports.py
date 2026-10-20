import importlib

from click.testing import CliRunner


def test_cli_help_smoke():
    cli_module = importlib.import_module("denovowest.cli")
    runner = CliRunner()

    result = runner.invoke(cli_module.main, ["--help"])

    assert result.exit_code == 0
    assert "simulation" in result.output
    assert "annotate-cadd" in result.output
    assert "annotate-custom" in result.output
    assert "annotate-dbnsfp" in result.output
    assert "annotate-vcf" in result.output
    assert "report" in result.output


def test_subcommand_help_smoke():
    cli_module = importlib.import_module("denovowest.cli")
    runner = CliRunner()

    for command_name in [
        "simulation",
        "annotate-cadd",
        "annotate-custom",
        "annotate-dbnsfp",
        "annotate-vcf",
        "report",
    ]:
        result = runner.invoke(cli_module.main, [command_name, "--help"])

        assert result.exit_code == 0, f"{command_name} help failed: {result.output}"
