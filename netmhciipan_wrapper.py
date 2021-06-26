import subprocess
import uuid
import sys


class NetMHCIIpanRun:
    def __init__(self, netmhc_path, allele, tmpdir_prefix, fasta):
        self.netmhc_path = netmhc_path
        self.allele = allele
        if tmpdir_prefix == "":
            tmpdir_prefix = "/tmp/"
        if not tmpdir_prefix.endswith("/"):
            tmpdir_prefix = tmpdir_prefix + "/"
        self.tmpdir = tmpdir_prefix + str(uuid.uuid4())
        self.fasta = fasta

    def run(self):
        args = [
            self.netmhc_path,
            "-filter",
            "-a", self.allele,
            "-tdir", self.tmpdir,
            "-f", self.fasta
        ]

        netmhc = subprocess.run(args, check=True, capture_output=True)
        if len(netmhc.stderr) > 0:
            raise subprocess.CalledProcessError(netmhc.stderr)
        results = self.parse_results(netmhc.stdout.decode(sys.stdout.encoding))
        return results

    def parse_results(self, netmhc_output):
        FILE_HEADER = 0
        HEADER = 1
        MIMICS = 2
        FOOTER = 3

        results = {}

        state = FILE_HEADER
        for line in netmhc_output.splitlines():
            if line.startswith("Error"):
                raise Exception("NetMHCIIPan Error:" + line)
            if state == "FILE_HEADER":
                if line == "":
                    continue
                elif line.startswith("#"):
                    continue
                else:
                    raise Exception("Bad Parsing: " + line)
            if line.startswith("-----"):
                state = (state + 1) % 4
                continue
            if state == HEADER:
                continue
            if state == MIMICS:
                ss = line.strip().split()
                pos, mhc, pep, of, core, core_rel, identity, score_el, rank_el, exp_bind, bind_level = ss
                if identity not in results:
                    results[identity] = {}
                if core not in results[identity]:
                    results[identity][core] = bind_level[2:]
                if core in results[identity]:
                    if results[identity][core] == "WB":
                        results[identity][core] = bind_level[2:]
            if state == FOOTER:
                continue

        return results
