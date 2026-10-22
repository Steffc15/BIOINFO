# Ex 3: send the top 10 codons to an LLM (Claude) and ask for a list
# of foods that contain the least of the associated amino acids.
# Uses the Claude Code CLI in headless mode (`claude -p`), so the request runs
# on the Claude subscription you are logged in with, not on paid API credits.
import glob
import os
import shutil
import subprocess

from ex0 import translate_codon
from ex2 import top_codons

AMINO_ACID_NAMES = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
    'C': 'Cysteine', 'E': 'Glutamic acid', 'Q': 'Glutamine', 'G': 'Glycine',
    'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
    'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
    'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
    '_': 'Stop codon',
}


def build_prompt(header, top10):
    lines = []
    for codon, count in top10:
        aa = translate_codon(codon)
        lines.append(f"- {codon} -> {AMINO_ACID_NAMES[aa]} ({aa}), found {count} times")
    amino_acids = sorted({AMINO_ACID_NAMES[translate_codon(c)] for c, _ in top10})
    return (
        f"I analysed the genome of {header} and found these 10 most frequent codons, "
        "each with its associated amino acid:\n"
        + "\n".join(lines)
        + "\n\nThe amino acids involved are: " + ", ".join(amino_acids) + ".\n"
        "Make a list of foods that contain the LOWEST amounts of these amino acids. "
        "For each amino acid give a few low-content foods, then a short combined list "
        "of foods that are low in all of them."
    )


def find_claude_cli():
    """Find the `claude` binary: on PATH, or the one shipped with the VS Code extension."""
    path = shutil.which("claude") or os.environ.get("CLAUDE_CODE_EXECPATH")
    if path:
        return path
    found = glob.glob(os.path.expanduser(
        "~/.vscode/extensions/anthropic.claude-code-*/resources/native-binary/claude"))
    if not found:
        raise FileNotFoundError("Claude Code CLI not found - install it or the VS Code extension.")
    return sorted(found)[-1]


def ask_claude(prompt):
    env = os.environ.copy()
    # Without an API key in the environment, the CLI uses the subscription login.
    env.pop("ANTHROPIC_API_KEY", None)
    result = subprocess.run([find_claude_cli(), "-p", prompt],
                            capture_output=True, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(result.stderr or result.stdout)
    return result.stdout


if __name__ == "__main__":
    header, _, top10 = top_codons()
    prompt = build_prompt(header, top10)
    print("PROMPT:\n" + prompt + "\n")

    print("LLM ANSWER:\n" + ask_claude(prompt))
