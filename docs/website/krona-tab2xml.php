<?php

$tabfile = "sunburst_allcats_filled.txt";
$lines = explode("\n", file_get_contents(rtrim($tabfile)));

$header = explode("\t", trim($lines[0]));

$levelcount = count($header) - 1;

$kxml = new SimpleXMLElement('<krona></krona>');

$nodes = array();

$attributes = $kxml -> addChild("attributes");
$attributes -> addAttribute("magnitude", "amount");
$att_amount = $attributes -> addChild("attribute", "amount");
$att_amount -> addAttribute("display", "amount");

$nodes["RootNode"] = $kxml -> addChild("node");
$nodes["RootNode"] -> addAttribute("name", "Leaf Protein");
$nodes["RootNode"] -> addChild("amount");
$nodes["RootNode"] -> amount = "0";

$cats = array();

//loop through lines in table and parse them to extract the nodeIDs and values
for ($r = 1; $r < count ($lines); $r++){

    if (substr_count($lines[$r], "\tend\t") > 0){

        $lineparts = explode("\tend\t", trim($lines[$r]));
        $val = floatval(array_pop(explode("\t", trim($lines[$r]))));

    } else if (substr_count($lines[$r], "\tNA\t") > 0){

        $lineparts = explode("\tNA\t", trim($lines[$r]));
        $val = floatval(array_pop(explode("\t", trim($lines[$r]))));

    }

    else {

        $line = explode("\t", trim($lines[$r]));
        $val = floatval(array_pop($line));
        $lineparts[0] = implode("\t", array_slice($line, 0, count($line)));

    }

    $catlevels = explode("\t", trim($lineparts[0]));

    //make a . delimited string ID for the node
    $catID = implode(".", $catlevels);

    if ($catID != "") { //just to make sure we don't have a blank line making an empty node
        for ($level = 0; $level < count($catlevels); $level++) {

            $node = array_slice($catlevels, 0, $level + 1);

            $nodeID = implode(".", $node);

            if (!array_key_exists($nodeID, $nodes)) {
//            echo ("<br>" . $nodeID);

                $nodename = $catlevels[$level];

                if ($level == 0) {

                    $nodes[$nodeID] = $nodes["RootNode"]->addChild("node");
                    $parentnodeID = "RootNode";
                    $nodes[$nodeID]->addAttribute("name", $nodename);
                    $nodes[$nodeID]->addAttribute("id", $nodeID);
                    $nodes[$nodeID]->addAttribute("level", $level);
                    $nodes[$nodeID]->addAttribute("parentID", $parentnodeID);
                    $nodes[$nodeID]->addChild("amount");
                    $nodes[$nodeID]->amount->addChild("val");
                    $nodes[$nodeID]->amount->val = $val;

                    $nodesbylevel[$level][$nodeID] = $nodes[$nodeID];

                } else {

                    $parentnode = array_slice($catlevels, 0, $level);
                    $parentnodeID = implode(".", $parentnode);

                    $nodes[$nodeID] = $nodes[$parentnodeID]->addChild("node");
                    $nodes[$nodeID]->addAttribute("name", $nodename);
                    $nodes[$nodeID]->addAttribute("id", $nodeID);
                    $nodes[$nodeID]->addAttribute("level", $level);
                    $nodes[$nodeID]->addAttribute("parentID", $parentnodeID);
                    $nodes[$nodeID]->addChild("amount");
                    $nodes[$nodeID]->amount->addChild("val");
                    $nodes[$nodeID]->amount->val = $val;

                    $nodesbylevel[$level][$nodeID] = $nodes[$nodeID];

                }
            }
        }
    }
}


//start at the deepest level and work up to the root
for ($level = $levelcount -1; $level > -1; $level--){

    //loop through all the nodes at the current level and add cumulatively add their value to their parent node's value
    foreach ($nodesbylevel[$level] as $nodeID => $node){

        $parentnodeID = $node['parentID'];
        $nodeval = floatval($nodes["$nodeID"]->amount->val);
        $parentnodeval = floatval($nodes["$parentnodeID"]->amount->val);
        $newparentnodeval = $nodeval + $parentnodeval;

        $nodes["$parentnodeID"]->amount->val = strval($newparentnodeval);

    }
}

//use XML string from SimpleXML object to create DomDocument object (which supports pretty XML formatting)
$xmldoc = new DomDocument();

$xmlstring = $kxml->asXML();
$xmldoc->loadXML($xmlstring);

$xmldoc->preserveWhiteSpace = false;
$xmldoc->formatOutput = true;

//generate pretty XML
$prettyxml = $xmldoc->saveXML();

//pretty XML is accessible by running script in browser then View Page Source
echo $prettyxml;


?>