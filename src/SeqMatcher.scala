/**
  * Created by Todd Stellanova on 20161026
  */

import java.io._
import scala.collection.mutable.Queue
import scala.io.Source
import util.control.Breaks._


object SeqMatcher {
      val filePath = "./data/coding_challenge_data_set.txt"
//    val filePath = "./data/simple_data.txt"
//    val filePath = "./data/remainder.txt"

  // KMP tables for input sentences, lazily created
  var kmpTables = scala.collection.mutable.Map[String, Array[Int]]()

  /**
    * Parse lines from a FASTA file
    * @param lines Iterator for lines in the FASTA file
    * @return A list of (name, sequence) from the lines
    */
  def parseFastaLines(lines: Iterator[String]): List[(String, String)] = {
    if (lines.isEmpty) return List()
    val name = lines.next.drop(1)
    val (seq, rest) = lines.span(_ (0) != '>')
    (name, seq.mkString) :: parseFastaLines(rest)
  }

  /**
    * Return a list of sentences (name, sequence)
    * from a FASTA formatted file
    * @param filePath
    * @return
    */
  def readFastaFile(filePath: String): List[(String, String)] = {
    parseFastaLines(Source.fromFile(filePath).getLines())
  }

  /**
    * Compute a Knuth–Morris–Pratt partial match table from a string
    * @see https://en.wikipedia.org/wiki/Knuth–Morris–Pratt_algorithm
    * @param s  The string from which to compute the KMP table
    * @return The KMP table
    */
  def computeKmpTable(s: String): Array[Int] = {
    val T = new Array[Int](s.length)
    var cnd = 0
    T(0) = -1
    T(1) = 0
    var pos = 2
    while (pos < s.length) {
      if (s(pos - 1) == s(cnd)) {
        T(pos) = cnd + 1
        pos += 1
        cnd += 1
      } else if (cnd > 0) {
        cnd = T(cnd)
      } else {
        T(pos) = 0
        pos += 1
      }
    }
    T
  }



  /*
  algorithm kmp_search:
    input:
  an array of characters, S (the text to be searched)
  an array of characters, W (the word sought)
  output:
    an integer (the zero-based position in S at which W is found)

  define variables:
    an integer, m ← 0 (the beginning of the current match in S)
  an integer, i ← 0 (the position of the current character in W)
  an array of integers, T (the table, computed elsewhere)

  while m + i < length(S) do
  if W[i] = S[m + i] then
  if i = length(W) - 1 then
  return m
  let i ← i + 1
  else
  if T[i] > -1 then
    let m ← m + i - T[i], i ← T[i]
  else
  let m ← m + 1, i ← 0

  (if we reach here, we have searched all of S unsuccessfully)
  return the length of S

  */

  /**
    * Search for needle in the base string, starting from 0 index
    * @param base The base string
    * @param baseMax The maximum base index to search
    * @param needle The string to search for
    * @param kmpTable  KMP partial match table
    * @return
    */
  def kmpSearch(base: String, baseMax: Int,
                needle: String, kmpTable: Array[Int]): (Int, Int) = {
    var m = 0 // the beginning of the current match in base
    var i = 0 // the position of the current character in needle

    while (( (m + i) < baseMax) && (i < needle.length)) {
      if (needle(i) == base(m + i)) {
        if (i == (needle.length - 1)) {
          return (m, i)
        }
        //count the matching chars
        i += 1
      }
      else {
        if (kmpTable(i) > -1) {
          m = m + i - kmpTable(i)
          i = kmpTable(i)
        }
        else {
          m = m + 1
          i = 0
        }
      }
    }
    if (i < 1) {
      m = -1
      i = 0
    }
    //if number characters matched > 0, return position and number
    (m, i)
  }


  /**
    * Append overlapping suffix to the base string at the given position,
    * ignoring suffixes that are already subsumed entirely by the base string
    * @param base
    * @param suffix
    * @param found Tuple of (start, length) where start is the position in the base string, length is suffix length
    * @return suffix appended to base string
    */
  def appendSuffix(base:String, suffix: String, found: (Int,Int) ): String = {
    val start = found._1
    val runLen = found._2
    val end = start+runLen

    if (end  < base.length) {
      //the base string is unmodified because
      //the found string falls entirely within the existing base
      println("subsumed " + end + " < " + base.length)
      return base
    }

    //it goes at the end
    base.substring(0, start) + suffix
  }


  /**
    * Try to find one overlapping sentence in the given Queue that we can glue to either the head or tail
    * of the given rootSequence
    * @param sentences Source sentences that have not yet been glued to the base sequence
    * @param remainders (out) Sentences that were tried by this method and no match was found
    * @param rootSequence The root sequence that we're trying to glue onto
    * @return Updated root sequence
    */
  def onePass(sentences: Queue[(String,String)], remainders: Queue[(String,String)], rootSequence:String):String  =  {
    var base = rootSequence
    //compute the KMP table (partial match table) for the root sequence
    val T = computeKmpTable(base)

    while (sentences.nonEmpty) {
      val curNode = sentences.dequeue()
      val curSeq = curNode._2
      var minSeqLen = curSeq.length
      if (base.length < curSeq.length)  minSeqLen = base.length

      //small optimization: Use smaller prefix of current sequence to find match
      var markerLen = minSeqLen / 48
      //This is tuned to example sentences provided but minimum could be adjusted
      if (markerLen < 2) markerLen = 2

      val curHead = curSeq.substring(0,markerLen)
      var detectLen = markerLen / 2
      if (detectLen < 1) detectLen = 1

      // is curSeq a suffix of base?
      val found = kmpSearch(base, base.length, curHead, T)
      if (found._2 >= detectLen) {
        //curSeq is a suffix of base: append curSeq to base
        println("suffix found  " + curNode._1 + " : " + found)
        base = appendSuffix(base,curSeq, (found._1, curSeq.length) )
        return base
      }
      else {
        if (!kmpTables.contains(curNode._1)) {
          //we've never computed the KMP table for this sequence before:
          //compute and cache it in case we need to reuse it multiple times
          kmpTables(curNode._1) = computeKmpTable(curSeq)
        }
        //is curSeq a prefix of base?
        val T2 = kmpTables(curNode._1)
        val baseHead = base.substring(0,markerLen)
        val found2 = kmpSearch(curSeq,curSeq.length, baseHead,T2)
        if (found2._2 > detectLen) {
          //curSeq is a prefix of base: prepend curSeq to base
          println("prefix found  " + curNode._1 + " : " + found2)
          base = curSeq.substring(0,found2._1) + base
          return base
        }
        else {
          //push to back of queue for another try later
          remainders.enqueue(curNode)
        }
      }
    }

    base
  }

  /**
    * Dump the queue of (name,sequence) as a "flattened" FASTA file,
    * where each sentence is on a single line
    * @param q
    * @param filePath
    */
  def dumpAsFlattenedFasta(q:Queue[(String, String)], filePath: String) = {
      val pw = new PrintWriter(new File(filePath))
      while (q.nonEmpty) {
        val curNode = q.dequeue()
        pw.println(">" + curNode._1)
        pw.println(curNode._2)
      }
      pw.close
  }

  def main(args: Array[String]) {
    println("Begin")

    var q1 = Queue[(String, String)]()
    var q2 = Queue[(String, String)]()
    //read all the sequences from file.
    //if we're concerned about RAM we could lazy load these
    println("readFastaFile from: " + filePath)
    var rawSequences = readFastaFile(filePath)
    println("num sentences: " + rawSequences.length)

    //stuff all the sentences into a queue
    q1 ++= rawSequences

//    dumpAsFlattenedFasta(q1, "./data/flattened.txt")
//    return

    //the root sequence could be head or tail or middle...doesn't matter
    var startNode = q1.dequeue()
    println("startNode: " + startNode._1)
    var rootSequence = startNode._2

    // the problem given states that overlapping strings overlap by >= half their length

    var lastCount = q1.length
    var curCount = -1
    breakable {
      while (q1.nonEmpty) {
        lastCount = curCount
        rootSequence = onePass(q1, q2, rootSequence)
        //      println("unglued: " + q2.length)
        q1 ++= q2
        println("unglued: " + q1.length)
        q2.clear()
        curCount = q1.length

        if (curCount == lastCount) {
          //we're stuck, unable to find any way to glue onto the rootSequence
          break
        }
      }
    }


    if (curCount > 0) {
      //we've failed to glue all the sentences onto the root sequence:
      //dump our remaining data
      val outFile = "./data/remainder.txt"
      println("FAILED! Writing remainders to " + outFile)
      q1.enqueue(("Root",rootSequence))
      dumpAsFlattenedFasta(q1, outFile)
    }
    else {
      println("SUCCESS!")
      println("result len: " + rootSequence.length)
      println(rootSequence)
      val pw = new PrintWriter(new File("./data/result.txt"))
      pw.println(rootSequence)
      pw.close
    }


  }

}
