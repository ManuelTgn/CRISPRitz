from ..crispritz_inputargs import CrispritzIndexingInputArgs
from .tst_manager import build_ternary_search_tree


def index_genome_cli(args: CrispritzIndexingInputArgs) -> None:
    build_ternary_search_tree(
        args.fastas, args.pam_file, args.bmax, args.outdir, args.threads, args.verbosity, args.debug
    )
