from ..crispritz_inputargs import CrispritzSearchInputArgs
from .tst_explorer import search_offtargets_tst

def search_offtargets_cli(args: CrispritzSearchInputArgs) -> None:
    search_offtargets_tst(args.indexes, args.pam_file, args.guides_file, args.mm, args.bdna, args.brna, args.outdir, args.threads, args.verbosity, args.debug)