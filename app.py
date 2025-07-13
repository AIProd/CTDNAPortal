import os, requests, streamlit as st, smtplib
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from typing import List, Dict
from Bio import Entrez
from langchain_openai import AzureChatOpenAI

# ── ENVIRONMENT VARIABLES ─────────────────────────────────────
ENTREZ_EMAIL      = "you@example.com"
AZURE_ENDPOINT    = os.environ["AZURE_OPENAI_ENDPOINT"]
AZURE_KEY         = os.environ["AZURE_OPENAI_KEY"]
BREVO_SENDER      = {"email": os.environ["BREVO_SENDER_EMAIL"],
                     "name":  os.environ["BREVO_SENDER_NAME"]}
BREVO_SMTP_PASS   = os.environ["BREVO_SMTP_PASSWORD"]

MAX_PAPERS_DEFAULT = 8
MAX_TRIALS_DEFAULT = 5
Entrez.email = ENTREZ_EMAIL

llm = AzureChatOpenAI(
    azure_endpoint     = AZURE_ENDPOINT,
    api_key            = AZURE_KEY,
    azure_deployment   = "gpt-4o",
    openai_api_version = "2024-12-01-preview",
    temperature        = 0,
)

# ── DATA FETCH FUNCTIONS ───────────────────────────────────────
def fetch_pubmed(term: str, max_results: int, lookback_days: int | None) -> List[Dict]:
    search_args = dict(db="pubmed", term=term, retmax=max_results)
    if lookback_days:
        since = (datetime.utcnow() - timedelta(days=lookback_days)).strftime("%Y/%m/%d")
        search_args.update(dict(mindate=since, datetype="pdat"))
    ids = Entrez.read(Entrez.esearch(**search_args))["IdList"]
    if not ids:
        return []
    docs = Entrez.read(Entrez.esummary(db="pubmed", id=",".join(ids)))
    return [
        {
            "title": doc.get("Title", "")[:250],
            "date":  doc.get("PubDate", "N/A"),
            "url":   f"https://pubmed.ncbi.nlm.nih.gov/{doc.get('Id')}/",
        }
        for doc in docs
    ]

def fetch_trials(condition: str, keyword: str, max_results: int) -> List[Dict]:
    r = requests.get(
        "https://clinicaltrials.gov/api/v2/studies",
        params={
            "query.cond": condition,
            "query.term": keyword,
            "pageSize":   max_results,
            "format":     "json",
        },
        timeout=20,
    )
    r.raise_for_status()
    studies = r.json().get("studies", [])
    results = []
    for s in studies:
        ps = s["protocolSection"]
        ident  = ps["identificationModule"]
        status = ps["statusModule"]
        phase  = ps.get("designModule", {}).get("phaseList", {}).get("phase", ["N/A"])[0]
        results.append(
            {
                "title": ident.get("officialTitle") or ident.get("briefTitle"),
                "nct":   ident["nctId"],
                "phase": phase,
                "start": status.get("startDateStruct", {}).get("startDate", "N/A"),
                "url":   f"https://clinicaltrials.gov/study/{ident['nctId']}",
            }
        )
    return results

def gpt_bullet(text: str) -> str:
    prompt = (
        "Summarize the following item for an oncology digest in ≤30 words, "
        "focusing on clinical relevance.\n\nItem:\n"
        f"{text}\n\nBullet:"
    )
    return llm.invoke(prompt).content.strip()

# ── EMAIL SENDER (via Brevo SMTP) ──────────────────────────────
def send_with_brevo(to_email: str, subject: str, body: str):
    msg = MIMEText(body)
    msg["Subject"] = subject
    msg["From"] = f"{BREVO_SENDER['name']} <{BREVO_SENDER['email']}>"
    msg["To"] = to_email

    with smtplib.SMTP("smtp-relay.brevo.com", 587) as server:
        server.starttls()
        server.login(BREVO_SENDER["email"], BREVO_SMTP_PASS)
        server.send_message(msg)

# ── STREAMLIT UI ───────────────────────────────────────────────
st.title("🧬 Oncology Digest Generator (via Brevo SMTP)")

topic = st.text_input("Search term", "ctDNA bladder cancer")
recipient = st.text_input("Recipient email")
max_papers = st.number_input("Max papers", 1, 100, MAX_PAPERS_DEFAULT, step=1)
max_trials = st.number_input("Max trials", 1, 50, MAX_TRIALS_DEFAULT, step=1)
lookback = st.selectbox("PubMed window",
                        ("Past 7 days", "Past 30 days", "All time"))
days_map = {"Past 7 days": 7, "Past 30 days": 30, "All time": None}

if st.button("Generate digest"):
    with st.spinner("Collecting literature…"):
        condition = "bladder cancer"
        keyword   = "ctDNA" if "ctDNA" in topic else topic

        papers = fetch_pubmed(topic, max_papers, days_map[lookback])
        trials = fetch_trials(condition, keyword, max_trials)

    for p in papers:
        p["summary"] = gpt_bullet(p["title"])
    for t in trials:
        t["summary"] = gpt_bullet(t["title"])

    today = datetime.now().strftime("%d %b %Y")
    lines = [
        f"Subject: Oncology Digest – {topic} ({today})",
        "",
        "Dear Dr. Bukavina,",
        "",
        f"Here is your digest on *{topic}*:",
        "",
        "📚 Papers:"
    ] + (
        sum([[f"• {p['summary']}", f"  ({p['date']}) {p['url']}"] for p in papers], [])
        or ["  – None found."]
    ) + [
        "",
        "🧪 Trials:"
    ] + (
        sum([[f"• {t['summary']} [{t['phase']}]",
              f"  NCT {t['nct']} — Start {t['start']}",
              f"  {t['url']}"] for t in trials], [])
        or ["  – None found."]
    ) + [
        "",
        "Regards,",
        "Oncology AI"
    ]
    digest_text = "\n".join(lines)

    st.text_area("Digest Preview", digest_text, height=400)

    if recipient:
        try:
            send_with_brevo(recipient, f"Oncology Digest – {topic}", digest_text)
            st.success(f"Email sent to {recipient}")
        except Exception as e:
            st.error(f"Email failed: {e}")
    else:
        st.warning("Please provide a recipient email to send.")
