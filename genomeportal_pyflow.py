#!/usr/bin/env python2.7
# coding=utf-8

"""
	Date:   12/11/2015
	Usage:  For usage instructions run with option --help
	Author: Madis Rumming <mrumming@cebitec.uni-bielefeld.de>
"""

__author__  = "Madis Rumming <mrumming@cebitec.uni-bielefeld.de>"
__copyright__ = "Copyright 2016, Computational Metagenomics, Faculty of Technology, Bielefeld University"

__version__ = "1.2.a"
__maintainer__ = "Madis Rumming"
__email__ = "mrumming@cebitec.uni-bielefeld.de"
__status__ = "Production"




import requests
import argparse
import os.path
import sys
from lxml import etree
import subprocess
import shlex

# Manual import of magic as module (DO NOT USE THE CEBITEC NATIVE ONE!!!)
import imp

magic = imp.load_source("magic", "/vol/cmg/share/virtualenvironments/pyflows/lib/python2.7/site-packages/magic.py")

# Extending PYTHONPATH for pyflow
sys.path.append("/vol/cmg/share/virtualenvironments/pyflows/lib/python2.7/site-packages/pyflow/")
from pyflow import WorkflowRunner


def parse_arguments():
    parser = argparse.ArgumentParser("Downloads automatically (meta)genomes from JGI genome portal in parallel.")
    parser.add_argument("-l", "--login-credentials", dest='loginfile',
                        help="Path to file containing JGI SSO login name in first line and password in second one",
                        required=True)
    parser.add_argument("-i", "--inputids", dest='inputids', help="IMG taxon OIDs as genome cart", required=True)
    parser.add_argument("-p", "--project-field", dest='project_field',
                        help="Index of project id in input file with oids, starting at 0.", required=False, default=-1,
                        type=int)
    parser.add_argument("--no-header", dest='skip_header', help="Set, if no header line is contained in genome cart",
                        required=False, default=True, action='store_false')

    parser.add_argument("--finished", dest='finished', help="Genome cart like file with finished IMG taxon OIDs",
                        required=False)
    parser.add_argument("-t", "--tmp-dir", dest='tmp_dir', help="TMP dir to use. Default: %s" % os.path.abspath('.'),
                        default=None, required=False)
    parser.add_argument("-f", "--final-dir", dest='dest_dir',
                        help="Base directory to save downloaded data (XML and downloaded files if desired) to. Default: %s" % os.path.abspath(
                            '.'), default=os.path.abspath('.'), required=False)

    parser.add_argument("--is-continued", dest='continued',
                        help="Enables continuing an erroneous or paused workflow. MUST use the same dataDirRoot as before.",
                        default=False, required=False, action='store_true')
    parser.add_argument("--is-dry-run", dest='dry_run', help="Check workflow without execution.", default=False,
                        action='store_true', required=False)

    parser.add_argument("-d", "--download-bundled", dest='download', help="Download bundled data", default=False,
                        action='store_true', required=False)
    parser.add_argument("-c", "--connection-limit", dest='con_limit', help="Connection limit. Default: %i" % (5),
                        type=int, default=5, required=False)

    parser.add_argument("-x", "--xml-dir", dest='xml_dir',
                        help="Directory to already downloaded XML descriptions. Use this, if you want to avoid being kicked by the IMG systems and already downloaded XML descriptions on (meta)genomes!",
                        default=None, type=str, required=False)

    parser.add_argument("-e", "--exclude-faa-fna-gff", dest='nofaafnagff',
                        help="Exclude directly fna, faa and gff files from .tar.gz", default=False, action='store_true',
                        required=False)

    parser.add_argument("-u", "--remove-unassembled", dest="unassembled",
                        help="Clean ALL files related to unassembled reads.", required=False, default=False,
                        action='store_true')
    parser.add_argument("-k", "--keep-unassembled-if-no-assembled", dest="keep_unassembled",
                        help="If -u set and only unassembled data is found (no assembled data is contained), these files are NOT removed (but fna, faa and gff are)!",
                        required=False, default=False, action='store_true')
    # parser.add_argument("--new-cluster", dest='new_cluster', help="Use the new cluster engine (OGE)", required=False, default=False, action='store_true')

    args = parser.parse_args()

    return (args)


class GatherXMLWorkflow(WorkflowRunner):
    def __init__(self, cur_ids, project_field, dest_dir, cookies):
        self.cur_ids = cur_ids
        self.project_field = project_field
        self.dest_dir = dest_dir
        self.cookies = cookies

    def workflow(self):
        for i in self.cur_ids:
            r = requests.get('http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=IMG_%s' % (i[0]),
                             cookies=self.cookies)
            if r.headers["Content-Type"] == "application/xml":
                outFile = open(os.path.join(self.dest_dir, "XML", "%s.xml" % (i[0])), "w")
                outFile.write(r.content)
                outFile.close()
            else:
                if not self.project_field == -1:
                    # GO VIA PROJECT ID -- http://genome.jgi.doe.gov/lookup?keyName=jgiProjectId&keyValue=407984
                    r = requests.get('http://genome.jgi.doe.gov/lookup?keyName=jgiProjectId&keyValue=%s' % (i[1]),
                                     cookies=self.cookies)
                    url = r.url.split("=")[-1]
                    r = requests.get('http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=%s' % (url),
                                     cookies=self.cookies)
                    if r.headers["Content-Type"] == "application/xml":
                        outFile = open(os.path.join(self.dest_dir, "XML", "proj_%s.xml" % (i[0])), "w")
                        outFile.write(r.content)
                        outFile.close()
                    else:
                        outFile = open(os.path.join(self.dest_dir, "XML", "ERR_%s.err" % (i[0])), "w")
                        outFile.write(r.content)
                        outFile.close()
                else:
                    outFile = open(os.path.join(self.dest_dir, "XML", "ERR_%s.err" % (i[0])), "w")
                    outFile.write(r.content)
                    outFile.close()


class GatherDownload(WorkflowRunner):
    def __init__(self, ids, cookies, dest_dir, tmp_dir, is_oid, omit, unassembled, keep_unassembled):
        self.ids = ids
        self.cookies = cookies
        self.dest_dir = dest_dir
        self.tmp_dir = tmp_dir
        self.is_oid = is_oid
        self.omit = omit
        self.unassembled = unassembled
        self.keep_unassembled = keep_unassembled

    def workflow(self):
        self.dest_dir = os.path.join(self.dest_dir, "Downloads")
        if not self.tmp_dir:
            self.tmp_dir = self.dest_dir

        if self.is_oid:
            for i in self.ids:
                r = requests.get('http://genome.jgi.doe.gov/IMG_%s/download/download_bundle.tar.gz' % (i),
                                 cookies=self.cookies)
                if not r.headers["Content-Type"] == "application/x-gzip":
                    outFile = open(os.path.join(self.dest_dir, "ERROR_%s.html" % (i)), "w")
                    outFile.write(r.content)
                    outFile.close()
                else:
                    if self.omit:
                        self.post_process_tar(r.content, i, self.tmp_dir, self.dest_dir, self.unassembled,
                                              self.keep_unassembled, "oid")
                    else:
                        print("WRITING FILE: %s" % (i))
                        outFile = open(os.path.join(self.dest_dir, "%s.tar.gz" % (i)), "w")
                        outFile.write(r.content)
                        outFile.close()

        else:
            for key in self.ids.keys():
                r = requests.get('http://genome.jgi.doe.gov/%s' % (self.ids[key]), cookies=self.cookies)
                if not r.headers["Content-Type"] == "application/x-gzip" and not r.headers[
                    "Content-Type"] == "application/octet-stream":

                    r = requests.get(
                        "https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=TaxonDetail&page=taxonDetail&taxon_oid=%s" % (
                            key), cookies=self.cookies)
                    root = etree.HTML(r.content)
                    refs = root.xpath(
                        ".//a[contains(@href, 'http://genome.jgi.doe.gov/lookup?keyName=jgiProjectId&keyValue=')]")
                    if len(refs) > 0:
                        proj_id = refs[0].attrib["href"].strip().split("=")[-1]

                        r = requests.get(
                            'http://genome.jgi.doe.gov/lookup?keyName=jgiProjectId&keyValue=%s' % (proj_id),
                            cookies=self.cookies)
                        url = r.url.split("=")[-1]
                        r = requests.get(
                            'http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=%s' % (url),
                            cookies=self.cookies)
                        if r.headers["Content-Type"] == "application/xml":
                            root = etree.XML(r.content)
                            files_ = root.findall(".//file[@url]")
                            url = ""
                            for f_ in files_:
                                if f_.attrib["filename"] == "download_bundle.tar.gz":
                                    url = f_.attrib["url"]
                            if not url:
                                for f_ in files_:
                                    if f_.attrib["filename"] == "%s.tar.gz" % (key):
                                        url = f_.attrib["url"]

                            if url:
                                r = requests.get(
                                    'http://genome.jgi.doe.gov/%s' % (url),
                                    cookies=self.cookies)
                                if not r.headers["Content-Type"] == "application/x-gzip" and not r.headers[
                                    "Content-Type"] == "application/octet-stream":
                                    outFile = open(os.path.join(self.dest_dir, "ERROR_%s.html" % (key)), "w")
                                    outFile.write(r.content)
                                    outFile.close()
                                else:
                                    if self.omit:
                                        self.post_process_tar(r.content, key, self.tmp_dir, self.dest_dir,
                                                              self.unassembled,
                                                              self.keep_unassembled, "proj2")
                                    else:
                                        print("WRITING FILE: %s" % (key))
                                        outFile = open(os.path.join(self.dest_dir, "%s.tar.gz" % (key)), "w")
                                        outFile.write(r.content)
                                        outFile.close()
                            else:
                                outFile = open(os.path.join(self.dest_dir, "ERROR_%s_proj.html" % (key)), "w")
                                outFile.write(r.content)
                                outFile.close()
                        else:
                            outFile = open(os.path.join(self.dest_dir, "ERROR_%s_proj.html" % (key)), "w")
                            outFile.write(r.content)
                            outFile.close()
                    else:
                        outFile = open(os.path.join(self.dest_dir, "ERROR_%s_proj.html" % (key)), "w")
                        outFile.write(r.content)
                        outFile.close()
                else:
                    if self.omit:
                        self.post_process_tar(r.content, key, self.tmp_dir, self.dest_dir, self.unassembled,
                                              self.keep_unassembled, "proj")
                    else:
                        print("WRITING FILE: %s" % (key))
                        outFile = open(os.path.join(self.dest_dir, "%s_proj.tar.gz" % (key)), "w")
                        outFile.write(r.content)
                        outFile.close()

    def post_process_tar(self, content, key, tmp_dir, dest_dir, unassembled, keep_unassembled, suffix):

        if suffix:
            suffix = "_%s" % (suffix)

        outFile = open(os.path.join(tmp_dir, "pre_%s%s.tar.gz" % (key, suffix)), "w+b")
        outFile.write(content)
        outFile.seek(0)

        outFinalFile = open(os.path.join(dest_dir, "shrinked", "%s%s.tar.gz" % (key, suffix)), "w")

        file_type = subprocess.check_output(['file', os.path.abspath(outFile.name)])
        file_type = file_type.split(" ", 1)[1].strip()

        try:

            cmd_tar_posix = "tar --delete --wildcards -f %s '*/*.f*a' '*/*.gff'"
            cmd_tar_gzip = "tar --delete --wildcards -f - '*/*.f*a' '*/*.gff'"

            if unassembled:
                if keep_unassembled:
                    cmd_tar_tf = shlex.split("tar tf %s" % (os.path.abspath(outFile.name)))
                    cont_files = subprocess.check_output(cmd_tar_tf)
                    cont_files = cont_files.splitlines()

                    for line in cont_files:
                        if ".a." in line or ".a," in line:
                            cmd_tar_posix += " '*/*.u*'"
                            cmd_tar_gzip += " '*/*.u*'"
                            break
                else:
                    cmd_tar_posix += " '*/*.u*'"
                    cmd_tar_gzip += " '*/*.u*'"

            if file_type.startswith("POSIX"):
                # gzipped tar



                cmd_tar = shlex.split(cmd_tar_posix % (os.path.abspath(outFile.name)))
                cmd_pigz_c = shlex.split("pigz -9 -p 24")

                tar_ = subprocess.Popen(cmd_tar, stdout=subprocess.PIPE)
                pigz_c = subprocess.check_call(cmd_pigz_c, stdin=tar_.stdout, stdout=outFinalFile)

                print(pigz_c)

            elif file_type.startswith("gzip"):
                # just tar


                cmd_pigz_d = shlex.split("pigz -d")
                cmd_tar = shlex.split(cmd_tar_gzip)
                cmd_pigz_c = shlex.split("pigz -9 -p 24")

                pigz_d = subprocess.Popen(cmd_pigz_d, stdin=outFile, stdout=subprocess.PIPE)
                tar_ = subprocess.Popen(cmd_tar, stdin=pigz_d.stdout, stdout=subprocess.PIPE)
                pigz_c = subprocess.check_call(cmd_pigz_c, stdin=tar_.stdout, stdout=outFinalFile)

                print(pigz_c)


            else:
                print("ERROR WITHIN %s" % (key))


        except:
            outFinalFile.close()
            outFile.close()
            cmd_except = shlex.split("rm %s" % (outFinalFile.name))
            subprocess.check_call(cmd_except)
            cmd_except = shlex.split(
                "mv %s %s" % (os.path.abspath(outFile.name), os.path.join(dest_dir, "%s_proj.tar.gz" % (key))))
            subprocess.check_call(cmd_except)
        finally:
            cmd_cleanup = shlex.split("rm %s" % (outFile.name))
            subprocess.check_call(cmd_cleanup)

            outFile.close()
            outFinalFile.close()


class GenomeportalWorkflow(WorkflowRunner):
    def __init__(self, oids, project_field, download_data, dest_dir, tmp_dir, login, pw, con_limit, xml_dir, omit,
                 unassembled, keep_unassembled):
        self.oids = oids
        self.project_field = project_field
        self.download_data = download_data
        self.dest_dir = dest_dir
        self.tmp_dir = tmp_dir
        self.login = login
        self.pw = pw
        self.con_limit = con_limit
        self.xml_dir = xml_dir
        self.omit = omit
        self.unassembled = unassembled
        self.keep_unassembled = keep_unassembled

    def workflow(self):

        if not os.path.exists(self.dest_dir):
            cmd = "mkdir -p %s" % (self.dest_dir)
            self.addTask(label="makeBaseDirectory", command=cmd, isForceLocal=True)
        else:
            self.addTask(label="makeBaseDirectory", isForceLocal=True)

        if not os.path.exists(os.path.join(self.dest_dir, "XML")):
            cmd = "mkdir %s" % (os.path.join(self.dest_dir, "XML"))
            self.addTask(label="makeXMLDirectory", command=cmd, isForceLocal=True, dependencies="makeBaseDirectory")
        else:
            self.addTask(label="makeXMLDirectory", isForceLocal=True)

        if not os.path.exists(os.path.join(self.dest_dir, "Downloads")):
            cmd = "mkdir %s" % (os.path.join(self.dest_dir, "Downloads"))
            self.addTask(label="makeDLDirectory", command=cmd, isForceLocal=True, dependencies="makeBaseDirectory")
        else:
            self.addTask(label="makeDLDirectory", isForceLocal=True)

        if self.xml_dir:
            if not os.path.exists(os.path.abspath(self.xml_dir)):
                sys.exit("Specified XML directory does not exist.")

        if self.tmp_dir:
            self.tmp_dir = os.path.join(os.path.abspath(self.tmp_dir), "gptdl")
            if not os.path.exists(self.tmp_dir):
                cmd = "mkdir -p %s" % (self.tmp_dir)
                self.addTask(label="makeTMPDirectory", command=cmd, isForceLocal=True)
            else:
                self.addTask(label="makeTMPDirectory", isForceLocal=True)
        else:
            self.addTask(label="makeTMPDirectory", isForceLocal=True)

        if self.omit:
            if not os.path.exists(os.path.join(os.path.abspath(self.dest_dir), "Downloads", "shrinked")):
                cmd = "mkdir %s" % (os.path.join(self.dest_dir, "Downloads", "shrinked"))
                self.addTask(label="makeDLDirectoryShrinked", command=cmd, isForceLocal=True,
                             dependencies="makeDLDirectory")
            else:
                self.addTask(label="makeDLDirectoryShrinked", isForceLocal=True)
        else:
            self.addTask(label="makeDLDirectoryShrinked", isForceLocal=True)

        os.popen(
            "curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=%s' --data-urlencode 'password=%s' -c cookies > /dev/null" % (
                self.login, self.pw))
        # cmd = "curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=%s' --data-urlencode 'password=%s' -c cookies > /dev/null" % (self.login, self.pw)
        # self.addTask(label="createCookie", command=cmd, isForceLocal=True)

        # self.waitForTasks("createCookie")

        inFile = open("cookies")
        cookies = {}

        for line in inFile:
            line = line.strip().split("\t")
            if line[0].startswith(".jgi"):
                cookies[line[-2]] = line[-1]
                break
        inFile.close()

        tasklist = []
        if not self.xml_dir:
            for i in xrange(self.con_limit):
                taskId = "XML%i" % (i)
                tasklist.append(taskId)
                cur_ids = self.oids[i::self.con_limit]
                print("Task:%s\t#IDs:%i" % (taskId, len(cur_ids)))
                wflow = GatherXMLWorkflow(cur_ids, self.project_field, self.dest_dir, cookies)
                self.addWorkflowTask(taskId, wflow,
                                     dependencies=["makeXMLDirectory", "makeDLDirectory", "makeTMPDirectory",
                                                   "makeDLDirectoryShrinked"])

        xml_dir = ""

        if not self.xml_dir:
            xml_dir = os.path.join(self.dest_dir, "XML")
        else:
            xml_dir = self.xml_dir

        tasklist2 = []
        self.waitForTasks(tasklist)
        if self.download_data:
            xmls = os.listdir(xml_dir)
            xmls = filter(lambda x: x.endswith(".xml"), xmls)
            xmls_oid = filter(lambda x: not x.startswith("proj"), xmls)
            xmls_oid = [x.split(".")[0] for x in xmls_oid]

            tasks_oids = []
            xmls_oid = [x[0] for x in self.oids if x[0] in xmls_oid]

            urls_oid = {}

            for oid in xmls_oid:
                if magic.from_file(os.path.join(xml_dir, "%s.xml" % (oid)), mime=True) == "application/xml":
                    print("Parsing file %s.xml" % (oid))
                    root = etree.parse(os.path.join(xml_dir, "%s.xml" % (oid)))
                    files_ = root.findall(".//file[@url]")
                    url = []
                    bundle = ""
                    for f_ in files_:
                        if f_.attrib["filename"] == "download_bundle.tar.gz":
                            bundle = f_.attrib["url"]
                    if bundle:
                        url.append(bundle)
                    else:
                        for f_ in files_:
                            if f_.attrib["filename"] == "%s.tar.gz" % (oid):
                                url.append(f_.attrib["url"])
                    if not len(url) == 0:
                        urls_oid[oid] = url[0]
                    else:
                        print("XML description, but no alternative .tar.gz available for %s" % (oid))
                        urls_oid[oid] = ""

            urls_keys = urls_oid.keys()
            for i in xrange(self.con_limit):
                urls_keys_ = urls_keys[i::self.con_limit]
                urls_ = {}
                for key in urls_keys_:
                    urls_[key] = urls_oid[key]
                taskId = "DLX%i" % (i)
                tasks_oids.append(taskId)
                wflow = GatherDownload(urls_, cookies, self.dest_dir, self.tmp_dir, False, self.omit,
                                       self.unassembled, self.keep_unassembled)
                self.addWorkflowTask(taskId, wflow, dependencies=tasklist)

            # for i in xrange(self.con_limit):
            #    taskId = "DL%i" % (i)
            #    tasks_oids.append(taskId)
            #    wflow = GatherDownload(xmls_oid[i::self.con_limit], cookies, self.dest_dir, self.tmp_dir, True,
            #                           self.omit, self.unassembled, self.keep_unassembled)
            #    self.addWorkflowTask(taskId, wflow, dependencies=tasklist)

            tasklist2.extend(tasks_oids)

            if not self.project_field == -1:
                xmls_proj = filter(lambda x: x.startswith("proj"), xmls)

                urls = {}

                for xml_path in xmls_proj:
                    if magic.from_file(os.path.join(xml_dir, xml_path), mime=True) == "application/xml":
                        print("Parsing file %s" % (xml_path))
                        root = etree.parse(os.path.join(xml_dir, xml_path))
                        files_ = root.findall(".//file[@url]")
                        url = []
                        xml_path = xml_path.split(".")[0].split("_")[1]
                        for f_ in files_:
                            if f_.attrib["filename"] == "%s.tar.gz" % (xml_path):
                                url.append(f_.attrib["url"])
                        if not len(url) == 0:
                            urls[xml_path] = url[0]
                        else:
                            print("XML description, but no alternative .tar.gz available for %s" % (xml_path))
                            urls[xml_path] = ""

                urls_keys = urls.keys()
                for i in xrange(self.con_limit):
                    urls_keys_ = urls_keys[i::self.con_limit]
                    urls_ = {}
                    for key in urls_keys_:
                        urls_[key] = urls[key]
                    taskId = "DLP%i" % (i)
                    tasklist2.append(taskId)
                    wflow = GatherDownload(urls_, cookies, self.dest_dir, self.tmp_dir, False, self.omit,
                                           self.unassembled, self.keep_unassembled)
                    self.addWorkflowTask(taskId, wflow, dependencies=tasks_oids)


                    ###
                    # NO CHECK AT ALL -> Will come in future release
                    # xmls_proj = [x.split(".")[0] for x in xmls_proj]
                    # xmls_proj = [x.split("_")[1] for x in xmls_proj]
                    # xmls_proj = [x for x in xmls_proj if x in self.oids]
                    #####

        cmd = "rm cookies"
        self.addTask(label="removeCookie", command=cmd, isForceLocal=True, dependencies=tasklist2)
        if self.tmp_dir:
            cmd = "rm -rf %s" % (self.tmp_dir)
            self.addTask(label="removeTMPDir", command=cmd, isForceLocal=True, dependencies="removeCookie")


def main():
    args = parse_arguments()

    ids = []

    # Export IMG taxonOIDs
    inFile = open(args.inputids)

    if args.skip_header:
        inFile.readline()

    for line in inFile:
        if not args.project_field == -1:
            l_ = []
            l_.append(line.strip().split('\t')[0])
            l_.append(line.strip().split('\t')[args.project_field])
            ids.append(l_)
        else:
            ids.append([line.strip().split('\t')[0]])
    inFile.close()

    login = ""
    pw = ""

    if not os.path.exists(os.path.abspath(args.loginfile)):
        sys.exit("File containing login credentials does not exist.")
    else:
        inFile = open(os.path.abspath(args.loginfile))
        login = inFile.readline().strip()
        pw = inFile.readline().strip()
        inFile.close()

    wflow = GenomeportalWorkflow(ids, args.project_field, args.download, args.dest_dir, args.tmp_dir, login, pw,
                                 args.con_limit, args.xml_dir, args.nofaafnagff, args.unassembled,
                                 args.keep_unassembled)
    retval = wflow.run(mode="sge", nCores="unlimited", memMb="unlimited", isDryRun=args.dry_run)
    sys.exit(retval)


if __name__ == "__main__":
    main()
